import requests
import time
from guerrillamail import GuerrillaMailSession
from bs4 import BeautifulSoup

from services import ss, batchtools

def get(seq):
	"""
	PROFsec - Part of PredictProtein suite
	Profile-based secondary structure prediction
	"""
	SS = ss.SS("PROFsec")
	SS.status = 0
	
	if (len(seq) < 10 or len(seq) > 4000):
		SS.pred += "Sequence length must be between 10 and 4000 amino acids"
		SS.conf += "Sequence length must be between 10 and 4000 amino acids"
		SS.status = 2
		print("PROFsec failed: Sequence length issue")
		return SS
	
	session = GuerrillaMailSession()
	email_address = session.get_session_state()['email_address']
	
	payload = {
		'sequence': seq,
		'email': email_address,
		'jobname': batchtools.randBase62(),
		'output': 'html'
	}
	
	try:
		# PredictProtein endpoint for PROFsec
		url = 'https://predictprotein.org/api/submit'
		r = requests.post(url, data=payload, timeout=30)
		
		if r.status_code != 200:
			SS.pred += f"Server returned HTTP {r.status_code}"
			SS.conf += f"Server returned HTTP {r.status_code}"
			SS.status = 2
			print(f"PROFsec failed: HTTP {r.status_code}")
			return SS
		
		soup = BeautifulSoup(r.text, 'html.parser')
		
		# Look for result URL or job ID
		result_link = None
		for link in soup.find_all('a', href=True):
			if 'result' in link.get('href', '').lower() or 'job' in link.get('href', '').lower():
				result_link = link['href']
				break
		
		if not result_link:
			# Try to parse results directly from response
			if 'PROF' in r.text or 'PROFsec' in r.text or 'prediction' in r.text.lower():
				lines = r.text.splitlines()
				for i, line in enumerate(lines):
					if ('PROF' in line or 'PROFsec' in line) and ('pred' in line.lower() or 'structure' in line.lower()):
						for j in range(i, min(i+10, len(lines))):
							if len(lines[j].strip()) == len(seq):
								SS.pred = lines[j].strip().replace('-', 'C')
								SS.status = 1
								print("PROFsec Complete")
								return SS
			
			SS.pred += "Could not parse results from server"
			SS.conf += "Could not parse results from server"
			SS.status = 2
			print("PROFsec failed: Could not parse results")
			return SS
		
		# Poll for results
		result_url = result_link if result_link.startswith('http') else url + result_link
		requesturl = batchtools.requestWait(result_url, 'PROFsec Not Ready', 20, 2700)
		
		if requesturl and requesturl.ok:
			raw = requesturl.text.splitlines()
			for i in range(len(raw)):
				if ('PROF' in raw[i] or 'PROFsec' in raw[i]) and len(raw[i].strip()) >= len(seq):
					SS.pred = raw[i].strip().replace('-', 'C')
					SS.status = 1
					print("PROFsec Complete")
					return SS
		
		SS.pred += "failed to respond after 45 minutes"
		SS.conf += "failed to respond after 45 minutes"
		SS.status = 2
		print("PROFsec failed: No response")
		
	except requests.RequestException as e:
		SS.pred += f"Network error: {str(e)}"
		SS.conf += f"Network error: {str(e)}"
		SS.status = 2
		print(f"PROFsec failed: {str(e)}")
	except Exception as e:
		SS.pred += f"Error: {str(e)}"
		SS.conf += f"Error: {str(e)}"
		SS.status = 2
		print(f"PROFsec failed: {str(e)}")
	
	return SS
