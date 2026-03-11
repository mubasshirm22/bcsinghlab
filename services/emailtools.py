import pickle
import email
import os
import os.path
import base64
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

SCOPES = ['https://www.googleapis.com/auth/gmail.modify']

#Logs in using token.pickle
def login():
	creds = None
	if os.path.exists('services/token.pickle'):
		with open('services/token.pickle', 'rb') as token:
			#print("opened pickle")
			creds = pickle.load(token)
		
	#Refresh if needed
	if creds and creds.expired and creds.refresh_token:
		creds.refresh(Request())

	service = build('gmail', 'v1', credentials=creds)
	
	return service

#Send email to target email. Email will be sent in HTML format, so normal text formatting will not work (such as \n)
def sendEmail(service, target, subject, msg):
	message = MIMEText(msg, 'html')
	message['to'] = target
	message['from'] = getEmailAddress(service)
	message['subject'] = subject
	message = {'raw': base64.urlsafe_b64encode(message.as_string().encode()).decode()}
	service.users().messages().send(userId='me', body=message).execute()

#Gets email address
def getEmailAddress(service):
	address = service.users().getProfile(userId='me').execute()
	return address['emailAddress']

#Returns the first email id that matches the query, -1 if not found
#Query guide: https://support.google.com/mail/answer/7190
#OR use the advanced search on the browser version of gmail and use that as the query
def searchEmailId(service, query):
	messageList = service.users().messages().list(userId='me',q=query).execute()
	
	if 'messages' in messageList:
		message = []
		message.extend(messageList['messages'])
		id = message[0]['id']
		message = service.users().messages().get(userId='me', id=id, format='raw').execute()
		return id
	return -1 #email not found

#Takes an email id and converts the 'raw' format of the email body to string
def decodeEmail(service, emailId):
	message = service.users().messages().get(userId='me', id=emailId, format='raw').execute()
	msg_str = base64.urlsafe_b64decode(message['raw'].encode('ASCII'))
	mime_msg = email.message_from_bytes(msg_str).as_string()
	return mime_msg
	
#Create token.pickle. Needed to change emails or update scopes
def createPickle():
	flow = InstalledAppFlow.from_client_secrets_file('services/credentials.json', SCOPES)
	creds = flow.run_local_server(port=0)
	# Save the credentials for the next run
	with open('services/token.pickle', 'wb') as token:
		pickle.dump(creds, token)


# ---------------------------------------------------------------------------
# Simple SMTP notification email — no OAuth required.
# Requires environment variables:
#   SMTP_EMAIL        — the Gmail address to send FROM (e.g. singhlab.notify@gmail.com)
#   SMTP_APP_PASSWORD — the 16-char Gmail App Password (not the account password)
# ---------------------------------------------------------------------------

def send_job_notification(to_email, job_id, site_url):
	"""
	Send a short notification email with a direct link to the results page.
	Returns True on success, False on failure (never raises).
	"""
	smtp_email = os.environ.get('SMTP_EMAIL')
	smtp_password = os.environ.get('SMTP_APP_PASSWORD')

	if not smtp_email or not smtp_password:
		print("Email notification skipped: SMTP_EMAIL / SMTP_APP_PASSWORD not set.")
		return False

	results_url = site_url.rstrip('/') + '/tools/sspred/dboutput/' + job_id

	html_body = """
	<div style="font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;max-width:520px;margin:0 auto;color:#0f172a;">
	  <div style="background:linear-gradient(135deg,#1e3a8a,#2563eb);padding:28px 32px;border-radius:10px 10px 0 0;">
	    <h2 style="color:#fff;margin:0;font-size:1.3rem;">SSPred — Prediction Complete</h2>
	    <p style="color:rgba(255,255,255,.8);margin:6px 0 0;font-size:0.9rem;">Singh Lab, Brooklyn College</p>
	  </div>
	  <div style="background:#fff;border:1px solid #e2e8f0;border-top:none;padding:28px 32px;border-radius:0 0 10px 10px;">
	    <p style="margin:0 0 16px;">Your secondary structure prediction job has finished.</p>
	    <p style="margin:0 0 8px;font-size:0.88rem;color:#64748b;">Job ID:</p>
	    <code style="font-family:monospace;background:#eff6ff;color:#1e3a8a;padding:4px 10px;border-radius:4px;font-size:0.9rem;">{job_id}</code>
	    <div style="margin:24px 0;">
	      <a href="{url}" style="display:inline-block;background:#2563eb;color:#fff;padding:12px 24px;border-radius:7px;text-decoration:none;font-weight:600;font-size:0.95rem;">View Results →</a>
	    </div>
	    <p style="font-size:0.82rem;color:#94a3b8;margin:0;">
	      You can also find this job in the <a href="{archive_url}" style="color:#2563eb;">Archive</a> using the Job ID above.<br>
	      This is an automated message from SSPred. Please do not reply.
	    </p>
	  </div>
	</div>
	""".format(job_id=job_id, url=results_url,
	           archive_url=site_url.rstrip('/') + '/tools/sspred/archive')

	try:
		msg = MIMEMultipart('alternative')
		msg['Subject'] = 'SSPred Job Complete — ' + job_id
		msg['From'] = 'SSPred <' + smtp_email + '>'
		msg['To'] = to_email
		msg.attach(MIMEText(
			'Your SSPred prediction job (' + job_id + ') is complete.\n\nView results: ' + results_url,
			'plain'
		))
		msg.attach(MIMEText(html_body, 'html'))

		with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
			server.login(smtp_email, smtp_password)
			server.sendmail(smtp_email, to_email, msg.as_string())

		print("Notification email sent to", to_email)
		return True
	except Exception as e:
		print("Failed to send notification email:", str(e))
		return False
