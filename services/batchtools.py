import time
import math
import requests
import secrets as _secrets
from guerrillamail import GuerrillaMailSession
import html
from bs4 import BeautifulSoup
#Contains functions related to output that are meant to be applied to multiple scripts

# ---------------------------------------------------------------------------
# Temporary email abstraction — tries GuerrillaMail first, falls back to
# mail.tm (free REST API, no account needed) if GuerrillaMail is down.
# ---------------------------------------------------------------------------

class _MailTmEmail:
	"""Minimal adapter so emailRequestWait can call .body on this object."""
	def __init__(self, body):
		self.body = body

class MailTmSession:
	"""
	Adapter that mimics the subset of GuerrillaMailSession used by the
	service files: session.get_email_list() and session.get_email(guid).body
	"""
	_MAILTM_BASE = 'https://api.mail.tm'

	def __init__(self, email_address, token):
		self.email_address = email_address
		self._token = token
		self._headers = {'Authorization': 'Bearer ' + token}

	def get_email_list(self):
		"""Returns list of objects with a .guid attribute (the message id)."""
		try:
			r = requests.get(self._MAILTM_BASE + '/messages', headers=self._headers, timeout=10)
			msgs = r.json().get('hydra:member', [])
			result = []
			for m in msgs:
				obj = type('_Msg', (), {'guid': m['id']})()
				result.append(obj)
			return result
		except Exception:
			return []

	def get_email(self, guid):
		"""Returns an object with a .body attribute (plain text body)."""
		try:
			r = requests.get(self._MAILTM_BASE + '/messages/' + guid, headers=self._headers, timeout=10)
			data = r.json()
			body = data.get('text', '') or data.get('html', '') or ''
			return _MailTmEmail(body)
		except Exception:
			return _MailTmEmail('')


def get_temp_email():
	"""
	Returns (email_address, session) using GuerrillaMail if available,
	mail.tm as fallback.  'session' is compatible with emailRequestWait.
	On total failure returns a static address with a None session
	(services that only need the address still work; inbox-reading services
	will time out gracefully).
	"""
	# Try GuerrillaMail first
	try:
		session = GuerrillaMailSession()
		email_address = session.get_session_state()['email_address']
		print("TempMail: using GuerrillaMail:", email_address)
		return email_address, session
	except Exception as e:
		print("TempMail: GuerrillaMail unavailable (" + str(e) + "), trying mail.tm ...")

	# Try mail.tm
	try:
		MAILTM = 'https://api.mail.tm'
		r = requests.get(MAILTM + '/domains', timeout=10)
		domain = r.json()['hydra:member'][0]['domain']
		username = 'sspred' + _secrets.token_hex(6)
		password = _secrets.token_hex(12)
		address = username + '@' + domain
		r2 = requests.post(MAILTM + '/accounts',
		                   json={'address': address, 'password': password},
		                   timeout=10)
		if r2.status_code in (200, 201):
			r3 = requests.post(MAILTM + '/token',
			                   json={'address': address, 'password': password},
			                   timeout=10)
			token = r3.json()['token']
			session = MailTmSession(address, token)
			print("TempMail: using mail.tm:", address)
			return address, session
	except Exception as e:
		print("TempMail: mail.tm unavailable (" + str(e) + "), using static fallback.")

	# Static fallback — services that only need to submit (PSI, JPred, Yaspin, etc.)
	# will still work fine. Inbox-reading services (Sable, SSPro, PSS) will time
	# out and report a clear failure message rather than crashing.
	static_addr = 'sspred.noreply@example.com'
	print("TempMail: using static fallback address:", static_addr)
	return static_addr, None

#Creates a random string to use for a prediction name. Can take a time and create a string from that
def randBase62(givenTime = None):
	if givenTime:
		integer = round(givenTime * 100000)
	else:
		integer = round(time.time() * 100000)
	chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
	result = ''
	while integer > 0:
		integer, remainder = divmod(integer, 62)
		result = chars[remainder]+result
	return result
	
#Takes current completed outputs and conducts a majority vote then returns it. If a majority vote results in an equal value, currently defaults to 'X'
def majorityVote(seq, ssObject):
	output = ''
	
	count = 0 #success counter
	for index in ssObject:
		if index.status == 1 or index.status == 3:
			count += 1
	
	
	if count >= 2: #vote only if at least 2 ssObjects are completed
		#create a counter for each character appearance
		seqLength = len(seq)
		cCount = [0] * seqLength
		hCount = [0] * seqLength
		eCount = [0] * seqLength
		
		for i in ssObject:
			if (i.status == 1 or i.status == 3) and len(i.pred) == seqLength:
				for j in range(0, seqLength):
					if i.pred[j] == 'C':
						cCount[j] += 1
					elif i.pred[j] == 'E':
						eCount[j] += 1
					elif i.pred[j] == 'H':
						hCount[j] += 1

		for i in range(0, seqLength):
			if eCount[i] > hCount[i] and eCount[i] > cCount[i]:
				output += 'E'
			elif hCount[i] > cCount[i] and hCount[i] > eCount[i]:
				output += 'H'
			elif cCount[i] > hCount[i] and cCount[i] > eCount[i]:
				output += 'C'
			else:
				output += 'X' #use X if unsure - typically shows when not all predictions are completed
	else:
		return None
	return output
import requests

def pdbget(pdbid, chain):
    # Validate input
    if not isinstance(pdbid, str) or not pdbid.isalnum():
        print("Invalid pdbid.")
        return None

    if not isinstance(chain, str) or not chain.isalnum():
        print("Invalid chain.")
        return None

    data = {
        "pdbid": pdbid,  # replace with pdbid variable
        "paste_field": "",
        "action": "compute",
        "contact_threshold": "6",
        "sensitive": "true"
    }
    
    url = 'https://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py'
    try:
        response = requests.post(url, data=data)
        response.raise_for_status()  # This will raise an exception for HTTP errors
    except requests.RequestException as e:
        print(f"Request failed: {e}")
        return None  # Ensure this is handled appropriately in your code

    lines = response.text.split('\n')

    selected_str_lines = []
    selected_seq_lines = []
    in_section_a = False
    for line in lines:
        if line.startswith('CHN') and f' {chain} ' in line:
            in_section_a = True
        elif line.startswith('CHN'):
            in_section_a = False
        if in_section_a:
            if line.startswith('STR'):
                selected_str_lines.append(line)
            elif line.startswith('SEQ'):
                selected_seq_lines.append(line)

    # If no sequence or structure data was found, return None
    if not selected_seq_lines or not selected_str_lines:
        print("No valid sequence or structure data found.")
        return None

    # Process sequences
    sequences = [line.split()[2] for line in selected_seq_lines if len(line.split()) > 2]
    sequence = ''.join(sequences)

    # Process structure, preserving the original logic for handling the structure information
    final_structure = ''.join([line[10:60] for line in selected_str_lines]).replace(" ", "C")

    # Construct the result dictionary
    result = {
        'pdbid': pdbid,
        'chain': chain,
        'primary': sequence,
        'color': [],  # Placeholder for color information
        'secondary': final_structure
    }

    return result


'''
#No auto canceling, infinite wait time

#Takes url to check, optional message for printing and optional sleep time in seconds. Defaults to 20 sec sleep time
#Returns the url when successful
def requestWait(requesturl, message = None, sleepTime = 20):
	while not requests.get(requesturl).ok:
		print(message)
		time.sleep(sleepTime)		
	return requests.get(requesturl)
	
#Takes a guerillamail session, search query, identifier line (Name: or Query:), and input name. Optional print message, and time to wait between checks
#Returns the bool email id and message when successful
def emailRequestWait(session, query, findLine, randName, printmsg = '', sleepTime = 60):
	message  = ''
	email_id = False
	
	while message == '': #loops until desired email is found or 15 min elapse
		print(printmsg)
		time.sleep(sleepTime)
		for e in session.get_email_list():			#For each email in inbox
			data = session.get_email(e.guid).body	#gets body of email
			if data is not None:					#Checks if email body is empty
				for dline in data.splitlines():		#Splits body into lines
					if findLine in dline:			#Checks if Query: line exists
						if dline[len(findLine):].strip() == randName:	#Checks if query is same as inputed seq name
							message = html.unescape(data)	#Sets message variable to email contents
							email_id = True			
	return email_id, message
'''

#Auto canceling versions
#Takes url to check, optional message for printing, and optional sleep time and cancel time in seconds. Defaults to 20 sec sleep time, 15 min wait to cancel
#Returns the url when successful
#Returns the url when successful
def requestWait(requesturl, message = None, sleepTime = 20 , cancelAfter = 1500):
	stime  = time.time()
	
	while not requests.get(requesturl).ok and time.time() < stime + cancelAfter: #loops until requesturl is found or cancelAfter min elapse
		print(message)
		time.sleep(sleepTime)
	return requests.get(requesturl)
	
#Takes a guerillamail session, search query, identifier line (Name: or Query:), and input name. Optional print message, time to wait between checks, and how long to wait until cancelling (both in seconds)
#Returns the bool email id and message when successful
def emailRequestWait(session, query, findLine, randName, printmsg = '', sleepTime = 15, cancelAfter = 1500):
	message  = ''
	stime = time.time()
	email_id = False
	
	while message == '' and time.time() < stime + cancelAfter: #loops until desired email is found or cancelAfter min elapse
		print(printmsg)
		time.sleep(sleepTime)
		try:
			print(session.get_email_list())
			for e in session.get_email_list():			#For each email in inbox
				data = session.get_email(e.guid).body	#gets body of email
				if data is not None:					#Checks if email body is empty
					for dline in data.splitlines():		#Splits body into lines
						if findLine in dline:			#Checks if Query: line exists
							if dline[len(findLine):].strip() == randName:	#Checks if query is same as inputed seq name
								message = html.unescape(data)	#Sets message variable to email contents
								email_id = True
		except:
			None
	return email_id, message