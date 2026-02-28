import requests
import time
from guerrillamail import GuerrillaMailSession

from services import ss, batchtools


def get(seq):
	print("DEBUG[PSI]: received sequence length:", len(seq))
	
	SS = ss.SS("PSI")

	if (len(seq) < 30 or len(seq) > 1500):
		if len(seq) < 30:
			SS.pred += "Sequence is too short (minimum 30 amino acids, got " + str(len(seq)) + ")"
			SS.conf += "Sequence is too short (minimum 30 amino acids, got " + str(len(seq)) + ")"
		else:
			SS.pred += "Sequence is too long (maximum 1500 amino acids, got " + str(len(seq)) + ")"
			SS.conf += "Sequence is too long (maximum 1500 amino acids, got " + str(len(seq)) + ")"
		SS.status = 2
		print("PsiPred failed: Sequence length issue - length:", len(seq))
		return SS 
		
	session = GuerrillaMailSession()
	email_address = session.get_session_state()['email_address']
	
	url = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submission/'

	# PSIPRED expects FASTA-formatted input. Wrap the raw sequence in a fixed
	# hidden FASTA header so the user can still just paste a plain sequence.
	fasta_seq = ">tr|B0R5N9|B0R5N9\n" + seq
	print("DEBUG[PSI]: FASTA format being sent (first 150 chars):", repr(fasta_seq[:150]))

	payload = {'input_data': fasta_seq}
	data = {
		'job': 'psipred',
		'submission_name': 'testing',
		'email': email_address,
	}
	
	try:
		r = requests.post(url, data=data, files=payload, headers={'accept': 'application/json'}, timeout=30)
		print("PSI submit response status:", r.status_code)
		print("PSI submit response headers:", dict(r.headers))
		try:
			print("PSI submit response body (truncated):", r.text[:500])
		except Exception:
			pass
		
		# Check HTTP status first
		if r.status_code != 200:
			error_msg = f"PSIPRED server returned HTTP {r.status_code}"
			if r.status_code == 500:
				error_msg += " (Internal Server Error - PSIPRED backend may be down)"
			elif r.status_code == 400:
				error_msg += " (Bad Request - check sequence format)"
			elif r.status_code == 503:
				error_msg += " (Service Unavailable)"
			SS.pred += error_msg
			SS.conf += error_msg
			SS.status = 2
			print("PsiPred failed:", error_msg)
			return SS
		
		# Try to parse JSON response
		try:
			resp_json = r.json()
		except ValueError as json_err:
			error_msg = f"PSIPRED returned invalid JSON. Response: {r.text[:200]}"
			SS.pred += error_msg
			SS.conf += error_msg
			SS.status = 2
			print("PsiPred failed:", error_msg, "JSON error:", str(json_err))
			return SS
		
		if 'UUID' not in resp_json:
			error_msg = f"PSIPRED response missing UUID. Full response: {str(resp_json)}"
			SS.pred += error_msg
			SS.conf += error_msg
			SS.status = 2
			print("PsiPred failed:", error_msg)
			return SS
		
		uuid = resp_json['UUID']
		print("PSI job submitted successfully, UUID:", uuid)

		jsonurl = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submission/' + uuid + '?format=json'
		r = requests.get(jsonurl, timeout=30)
		
		if r.status_code != 200:
			error_msg = f"Failed to check PSIPRED job status (HTTP {r.status_code})"
			SS.pred += error_msg
			SS.conf += error_msg
			SS.status = 2
			print("PsiPred failed:", error_msg)
			return SS

		filesUUID = r.json()['submissions'][0]['UUID'] 
		horiz = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submissions/' + filesUUID + '.horiz'
		
		#Length 1500 takes around 5 min, increased timeout to 45 min
		requesturl = batchtools.requestWait(horiz, 'PsiPred Not Ready', 20, 2700)
		
		if requesturl and requesturl.ok:
			raw = requesturl.text.splitlines()
			for i in range(len(raw)):
				raw[i] = raw[i].strip()
				if raw[i].startswith("Conf"):
					SS.conf += raw[i][6:]
				if raw[i].startswith("Pred"):
					SS.pred += raw[i][6:]
					
			SS.status = 1
			print("PsiPred Complete")
		else:
			SS.pred += "PSIPRED job timed out after 45 minutes (results not available)"
			SS.conf += "PSIPRED job timed out after 45 minutes (results not available)"
			SS.status = 2
			print("PsiPred failed: No response after 45 minutes")

	except requests.RequestException as req_err:
		error_msg = f"Network error contacting PSIPRED: {str(req_err)}"
		SS.pred += error_msg
		SS.conf += error_msg
		SS.status = 2
		print("PsiPred failed: Network error:", str(req_err))
	except Exception as e:
		error_msg = f"PSIPRED error: {type(e).__name__}: {str(e)}"
		SS.pred += error_msg
		SS.conf += error_msg
		SS.status = 4
		print("PsiPred failed:", error_msg)
		import traceback
		traceback.print_exc()
			
	print("PSI::")
	print(SS.pred)
	print(SS.conf)
	
	return SS
	