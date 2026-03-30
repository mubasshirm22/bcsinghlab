import io
import re
import time
import requests

from services import ss

_WEBFACE_URL  = 'https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi'
_CONFIG_FILE  = '/var/www/services/services/NetSurfP-2.0/webface.cf'
_POLL_SLEEP   = 20   # seconds between poll attempts
_CANCEL_AFTER = 1500 # 25 minutes total


def get(seq):
	SS = ss.SS("NetSurf")
	SS.status = 0

	if len(seq) < 10 or len(seq) > 4000:
		SS.pred   = "Sequence length must be between 10 and 4000 amino acids"
		SS.conf   = "Sequence length must be between 10 and 4000 amino acids"
		SS.status = 2
		print("NetSurf failed: sequence length out of range")
		return SS

	# ── 1. Submit job ──────────────────────────────────────────────────────
	fasta_bytes = (">query\n" + seq + "\n").encode()

	files = {
		'uploadfile': ('query.fasta', io.BytesIO(b''), 'application/octet-stream'),
	}
	data = {
		'configfile': _CONFIG_FILE,
		'pastefile':  ">query\n" + seq + "\n",
	}

	try:
		r = requests.post(
			_WEBFACE_URL,
			data=data,
			files=files,
			allow_redirects=True,
			timeout=60,
		)
	except requests.RequestException as e:
		SS.pred   = "Network error during submission: " + str(e)
		SS.conf   = SS.pred
		SS.status = 2
		print("NetSurf failed (submit):", e)
		return SS

	# Extract jobid from the final URL or response body
	jobid = _extract_jobid(r)
	if not jobid:
		SS.pred   = "Could not obtain job ID from NetSurf server"
		SS.conf   = SS.pred
		SS.status = 2
		print("NetSurf failed: no jobid found. Final URL:", r.url)
		print("NetSurf response (first 500):", r.text[:500])
		return SS

	print("NetSurf jobid:", jobid)

	# ── 2. Poll until done ─────────────────────────────────────────────────
	csv_text = _poll(jobid)
	if csv_text is None:
		SS.pred   = "NetSurf job timed out or returned no result"
		SS.conf   = SS.pred
		SS.status = 2
		print("NetSurf failed: poll timed out")
		return SS

	# ── 3. Parse CSV → pred string ─────────────────────────────────────────
	pred = _parse_csv(csv_text, len(seq))
	if pred is None:
		SS.pred   = "Could not parse NetSurf prediction output"
		SS.conf   = SS.pred
		SS.status = 2
		print("NetSurf failed: CSV parse error")
		print("NetSurf CSV (first 500):", csv_text[:500])
		return SS

	SS.pred   = pred
	SS.conf   = "No confidence (NetSurf)"
	SS.status = 3   # status 3 = valid pred but no numeric confidence
	print("NetSurf complete, pred length:", len(SS.pred))
	print("NETSURF::")
	print(SS.pred)
	return SS


# ── helpers ────────────────────────────────────────────────────────────────

def _extract_jobid(response):
	"""Pull jobid from redirect URL or response HTML."""
	# Try final URL query param  ?jobid=XXXX
	m = re.search(r'[?&]jobid=([A-Za-z0-9_\-]+)', response.url)
	if m:
		return m.group(1)
	# Try body
	m = re.search(r'[?&]jobid=([A-Za-z0-9_\-]+)', response.text)
	if m:
		return m.group(1)
	# Some servers embed it as  jobid = "XXXX"
	m = re.search(r'jobid\s*[=:]\s*["\']?([A-Za-z0-9_\-]{6,})["\']?', response.text, re.IGNORECASE)
	if m:
		return m.group(1)
	return None


def _poll(jobid):
	"""Poll the ajax endpoint until the result CSV is available.
	Returns the CSV text, or None on timeout."""
	deadline = time.time() + _CANCEL_AFTER
	while time.time() < deadline:
		time.sleep(_POLL_SLEEP)
		try:
			r = requests.get(
				_WEBFACE_URL,
				params={'ajax': '1', 'jobid': jobid, 'wait': '20'},
				timeout=60,
			)
		except requests.RequestException as e:
			print("NetSurf poll error:", e)
			continue

		text = r.text.strip()
		print("NetSurf poll response (first 200):", text[:200])

		# If it looks like CSV data with H/E/C column, we're done
		if _looks_like_csv(text):
			return text

		# webface2 sometimes returns a redirect URL to the actual result file
		if r.status_code in (301, 302, 303):
			loc = r.headers.get('Location', '')
			if loc:
				try:
					r2 = requests.get(loc, timeout=60)
					if _looks_like_csv(r2.text):
						return r2.text
				except Exception:
					pass

		# Still running — keep waiting
	return None


def _looks_like_csv(text):
	"""Heuristic: NetSurf CSV has lines with comma-separated fields,
	and the q3 column contains H, E, or C."""
	if not text or len(text) < 20:
		return False
	# Must have at least one non-comment line with H, E, or C in expected position
	for line in text.splitlines():
		if line.startswith('#') or not line.strip():
			continue
		parts = line.split(',')
		if len(parts) >= 4 and parts[3].strip() in ('H', 'E', 'C'):
			return True
	return False


def _parse_csv(csv_text, expected_len):
	"""Parse NetSurf-P 2.0 CSV into an H/E/C prediction string.

	CSV columns (0-indexed):
	  0=id, 1=n (1-based residue index), 2=seq, 3=q3, 4=q8, 5=rsa, ...
	"""
	residues = {}
	for line in csv_text.splitlines():
		line = line.strip()
		if not line or line.startswith('#'):
			continue
		parts = line.split(',')
		if len(parts) < 4:
			continue
		try:
			n   = int(parts[1].strip())
			q3  = parts[3].strip().upper()
		except (ValueError, IndexError):
			continue
		if q3 in ('H', 'E', 'C'):
			residues[n] = q3

	if not residues:
		return None

	# Build string from residue 1 … max_index
	max_idx = max(residues.keys())
	pred = ''.join(residues.get(i, 'C') for i in range(1, max_idx + 1))

	# If parsed length doesn't match expected, something went wrong
	if abs(len(pred) - expected_len) > 5:
		print(f"NetSurf length mismatch: expected {expected_len}, got {len(pred)}")
		return None

	return pred
