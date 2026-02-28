import requests
import time
from bs4 import BeautifulSoup
import re

from services import ss, batchtools

def get(seq):
	"""
	Predator - Secondary structure prediction using sequence alignment
	Server: npsa-prabi.ibcp.fr
	Endpoint: https://npsa-prabi.ibcp.fr/cgi-bin/secpred_preda.pl
	"""
	SS = ss.SS("Predator")
	SS.status = 0
	
	if (len(seq) < 10 or len(seq) > 4000):
		SS.pred += "Sequence length must be between 10 and 4000 amino acids"
		SS.conf += "Sequence length must be between 10 and 4000 amino acids"
		SS.status = 2
		print("Predator failed: Sequence length issue")
		return SS
	
	payload = {
		'title': '',  # Optional, can be empty
		'notice': seq,  # The protein sequence
		'ali_width': '70',  # Alignment width
		'predatorssmat': 'dssp'  # Secondary structure matrix type
	}
	
	try:
		# Correct Predator endpoint from browser inspection
		url = 'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_preda.pl'
		
		# POST request with form data (allow more time; Predator can be slow)
		r = requests.post(url, data=payload, timeout=180, headers={
			'Content-Type': 'application/x-www-form-urlencoded'
		})
		
		if r.status_code != 200:
			SS.pred += f"Server returned HTTP {r.status_code}"
			SS.conf += f"Server returned HTTP {r.status_code}"
			SS.status = 2
			print(f"Predator failed: HTTP {r.status_code}")
			return SS
		
		# Parse the HTML response
		soup = BeautifulSoup(r.text, 'html.parser')
		
		# Try multiple strategies to find the prediction
		prediction_found = None
		
		# Strategy 1: Look for <pre> tags (common for formatted output)
		for pre in soup.find_all('pre'):
			text = pre.get_text().strip()
			# Check if it's a valid secondary structure string
			if len(text) == len(seq) and all(c in 'HEC-' for c in text.upper()):
				prediction_found = text
				break
		
		# Strategy 2: Look in table cells
		if not prediction_found:
			for td in soup.find_all('td'):
				text = td.get_text().strip()
				if len(text) == len(seq) and all(c in 'HEC-' for c in text.upper()):
					prediction_found = text
					break
		
		# Strategy 3: Look for text that matches sequence length in the raw HTML
		if not prediction_found:
			# Get all text and split by lines
			text_content = soup.get_text()
			lines = text_content.splitlines()
			
			for line in lines:
				line_stripped = line.strip()
				# Check if line matches sequence length and contains only H/E/C/-
				if len(line_stripped) == len(seq) and all(c in 'HEC-' for c in line_stripped.upper()):
					prediction_found = line_stripped
					break
		
		# Strategy 4: Look for patterns like "Pred:" or "Structure:" followed by prediction
		if not prediction_found:
			text_content = soup.get_text()
			# Look for patterns like "Pred:" or "Structure:" followed by the prediction
			patterns = [
				r'Pred[iction]*:?\s*([HEC\-\s]+)',
				r'Structure:?\s*([HEC\-\s]+)',
				r'Secondary[_\s]?Structure:?\s*([HEC\-\s]+)',
			]
			for pattern in patterns:
				match = re.search(pattern, text_content, re.IGNORECASE)
				if match:
					pred_str = match.group(1).strip().replace(' ', '').replace('\n', '')
					if len(pred_str) == len(seq) and all(c in 'HEC-' for c in pred_str.upper()):
						prediction_found = pred_str
						break
		
		# Strategy 5: Check if results are in a specific div or section
		if not prediction_found:
			# Look for divs with class/id containing "result", "pred", "structure"
			for div in soup.find_all(['div', 'span', 'p']):
				class_attr = ' '.join(div.get('class', [])) if div.get('class') else ''
				id_attr = div.get('id', '')
				if any(keyword in (class_attr + ' ' + id_attr).lower() for keyword in ['result', 'pred', 'structure', 'output']):
					text = div.get_text().strip()
					# Try to extract prediction from this section
					for line in text.splitlines():
						line_stripped = line.strip()
						if len(line_stripped) == len(seq) and all(c in 'HEC-' for c in line_stripped.upper()):
							prediction_found = line_stripped
							break
					if prediction_found:
						break
		
		# If we found a prediction, use it
		if prediction_found:
			SS.pred = prediction_found.replace('-', 'C')
			SS.status = 1
			SS.conf = "Predator Does Not Provide Confidence"
			print("Predator Complete")
			return SS
		
		# Check if there's a message about processing/queuing
		text_content = soup.get_text().lower()
		if 'processing' in text_content or 'queued' in text_content or 'wait' in text_content:
			# Results might be processing - check for result URL to poll
			result_link = None
			for link in soup.find_all('a', href=True):
				href = link.get('href', '')
				if any(keyword in href.lower() for keyword in ['result', 'job', 'status', 'output']):
					result_link = href
					break
			
			if result_link:
				# Poll for results
				result_url = result_link if result_link.startswith('http') else 'https://npsa-prabi.ibcp.fr' + result_link
				requesturl = batchtools.requestWait(result_url, 'Predator Not Ready', 20, 2700)
				
				if requesturl and requesturl.ok:
					# Parse the final results with same strategies
					soup_result = BeautifulSoup(requesturl.text, 'html.parser')
					
					# Try all parsing strategies again on the result page
					for pre in soup_result.find_all('pre'):
						text = pre.get_text().strip()
						if len(text) == len(seq) and all(c in 'HEC-' for c in text.upper()):
							SS.pred = text.replace('-', 'C')
							SS.status = 1
							SS.conf = "Predator Does Not Provide Confidence"
							print("Predator Complete")
							return SS
					
					# Try table cells
					for td in soup_result.find_all('td'):
						text = td.get_text().strip()
						if len(text) == len(seq) and all(c in 'HEC-' for c in text.upper()):
							SS.pred = text.replace('-', 'C')
							SS.status = 1
							SS.conf = "Predator Does Not Provide Confidence"
							print("Predator Complete")
							return SS
					
					# Try lines
					text_content = soup_result.get_text()
					lines = text_content.splitlines()
					for line in lines:
						line_stripped = line.strip()
						if len(line_stripped) == len(seq) and all(c in 'HEC-' for c in line_stripped.upper()):
							SS.pred = line_stripped.replace('-', 'C')
							SS.status = 1
							SS.conf = "Predator Does Not Provide Confidence"
							print("Predator Complete")
							return SS
			
			SS.pred += "Predator server accepted the job but results could not be retrieved (format or URL may have changed)."
			SS.conf += "Predator server accepted the job but results could not be retrieved."
			SS.status = 2
			print("Predator failed: Could not retrieve processing results")
			return SS
		
		# If we get here, couldn't find prediction
		SS.pred += "Could not parse prediction from Predator server response (results may still be available on the website)."
		SS.conf += "Could not parse prediction from Predator server response."
		SS.status = 2
		print("Predator failed: Could not parse results")
		print("Response preview (first 500 chars):", r.text[:500])
		
	except requests.RequestException as e:
		SS.pred += f"Network error: {str(e)}"
		SS.conf += f"Network error: {str(e)}"
		SS.status = 2
		print(f"Predator failed: {str(e)}")
	except Exception as e:
		SS.pred += f"Error: {str(e)}"
		SS.conf += f"Error: {str(e)}"
		SS.status = 2
		print(f"Predator failed: {str(e)}")
		import traceback
		traceback.print_exc()
	
	return SS
