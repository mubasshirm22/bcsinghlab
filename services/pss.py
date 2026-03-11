import requests
from bs4 import BeautifulSoup
import time
from guerrillamail import GuerrillaMailSession

from services import ss, batchtools

def get(seq):

	SS = ss.SS("PSS")
	SS.status = 0

	if (len(seq) > 4000):
		SS.pred += "Sequence longer than 4000"
		SS.conf += "Sequence longer than 4000"
		SS.status = 2 #error status
		print("PSSPred failed: Sequence longer than 4000")
		return SS #return SS so it will be readable as an ssObject
		
	email_address, session = batchtools.get_temp_email()

	randName = batchtools.randBase62()
	fasta_seq = '>testprot\n' + seq

	payload = {'REPLY-E-MAIL': email_address,
		'TARGET-NAME': randName,
		'SEQUENCE': fasta_seq}
	files = {'seq_file': ('', b'', 'application/octet-stream')}
	r = requests.post('https://aideepmed.com/PSSpred/bin/receive.cgi', data=payload, files=files)

	print("PSSPred submit status:", r.status_code)
	print("PSSPred submit response:", r.text[:300])

	# Results are emailed; wait for email containing the job name
	query = 'from:(PSSpred) subject:(PSSpred) ' + randName
	email_id, message = batchtools.emailRequestWait(session, query, "Name: ", randName, "PSSpred Not Ready", 30, 2700)

	if email_id:
		raw = message.splitlines()
		for i in range(len(raw)):
			if raw[i].startswith("conf"):
				SS.conf += raw[i][6:].strip()
			if raw[i].startswith("SS"):
				SS.pred += raw[i][6:].strip()

		SS.status = 1
		print("PSSPred Complete")
	else:
		SS.pred += "failed to respond in time"
		SS.conf += "failed to respond in time"
		SS.status = 2 #error status
		print("PSSPred failed: No response")
	return SS