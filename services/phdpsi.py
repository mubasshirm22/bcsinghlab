import requests
import time
from guerrillamail import GuerrillaMailSession
from bs4 import BeautifulSoup

from services import ss, batchtools

def get(seq):
	"""
	PHDpsi - Part of PredictProtein suite
	Neural network-based secondary structure prediction
	"""
	SS = ss.SS("PHDpsi")
	SS.status = 0
	
	if (len(seq) < 10 or len(seq) > 4000):
		SS.pred += "Sequence length must be between 10 and 4000 amino acids"
		SS.conf += "Sequence length must be between 10 and 4000 amino acids"
		SS.status = 2
		print("PHDpsi failed: Sequence length issue")
		return SS
	
	email_address, _session = batchtools.get_temp_email()
	
	# PredictProtein typically uses this endpoint
	# Note: This may need adjustment based on actual server response
	payload = {
		'sequence': seq,
		'email': email_address,
		'jobname': batchtools.randBase62(),
		'output': 'html'
	}
	
	try:
		# PredictProtein main endpoint - may need to be adjusted
		url = 'https://predictprotein.org/api/submit'
		r = requests.post(url, data=payload, timeout=30)
		
		if r.status_code != 200:
			SS.pred += f"Server returned HTTP {r.status_code}"
			SS.conf += f"Server returned HTTP {r.status_code}"
			SS.status = 2
			print(f"PHDpsi failed: HTTP {r.status_code}")
			return SS
		
		# Try to parse response for job ID or results
		soup = BeautifulSoup(r.text, 'html.parser')
		
		# Look for result URL or job ID in response
		# This pattern may need adjustment based on actual server response
		result_link = None
		for link in soup.find_all('a', href=True):
			if 'result' in link.get('href', '').lower() or 'job' in link.get('href', '').lower():
				result_link = link['href']
				break
		
		if not result_link:
			# Try alternative: check if response contains results directly
			if 'PHD' in r.text or 'prediction' in r.text.lower():
				# Parse results from HTML
				lines = r.text.splitlines()
				for i, line in enumerate(lines):
					if 'PHD' in line and ('pred' in line.lower() or 'structure' in line.lower()):
						# Extract prediction - pattern may vary
						for j in range(i, min(i+10, len(lines))):
							if len(lines[j].strip()) == len(seq):
								SS.pred = lines[j].strip().replace('-', 'C')
								SS.status = 1
								print("PHDpsi Complete")
								return SS
			
			SS.pred += "Could not parse results from server"
			SS.conf += "Could not parse results from server"
			SS.status = 2
			print("PHDpsi failed: Could not parse results")
			return SS
		
		# Poll for results if we got a job ID
		result_url = result_link if result_link.startswith('http') else url + result_link
		requesturl = batchtools.requestWait(result_url, 'PHDpsi Not Ready', 20, 2700)
		
		if requesturl and requesturl.ok:
			raw = requesturl.text.splitlines()
			for i in range(len(raw)):
				if 'PHD' in raw[i] and len(raw[i].strip()) >= len(seq):
					SS.pred = raw[i].strip().replace('-', 'C')
					SS.status = 1
					print("PHDpsi Complete")
					return SS
		
		SS.pred += "failed to respond after 45 minutes"
		SS.conf += "failed to respond after 45 minutes"
		SS.status = 2
		print("PHDpsi failed: No response")
		
	except requests.RequestException as e:
		SS.pred += f"Network error: {str(e)}"
		SS.conf += f"Network error: {str(e)}"
		SS.status = 2
		print(f"PHDpsi failed: {str(e)}")
	except Exception as e:
		SS.pred += f"Error: {str(e)}"
		SS.conf += f"Error: {str(e)}"
		SS.status = 2
		print(f"PHDpsi failed: {str(e)}")
	
	return SS
