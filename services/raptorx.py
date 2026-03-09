import requests
import time
from bs4 import BeautifulSoup
from lxml import html
from zipfile import ZipFile
from urllib import request
import re
import os
import io
from guerrillamail import GuerrillaMailSession

from services import ss, batchtools


def get(seq):

	SS = ss.SS("RaptorX")
	SS.status = 0
	
	if (len(seq) < 27 or len(seq) > 4000): #<27 shouldnt happen with input validation
		SS.pred += "Sequence is longer than 4000"
		SS.conf += "Sequence is longer than 4000"
		SS.status = 2 #error status
		print("RaptorX failed: Sequence is shorter than 27 or longer than 4000")
		'''
		SS.pred += "Sequence is shorter than 27 or longer than 4000"
		SS.conf += "Sequence is shorter than 27 or longer than 4000"
		SS.status = 2 #error status
		print("RaptorX failed: Sequence is shorter than 27 or longer than 4000")
		'''
		return SS #return SS so it will be readable as an ssObject
	
	session = GuerrillaMailSession()	#Creates GuerrillaMail session
	email_address = session.get_session_state()['email_address'] #retrieves temp email address

	payload = {'jobname': 'myprot', 
		'useProfile': 'true', 
		'predict_sub': 'Submit', 
		'email': email_address, 
		'sequences': seq}

	r = requests.post('http://raptorx.uchicago.edu/StructurePropertyPred/predict/pred/', data=payload)


	soup = BeautifulSoup(r.text,'html.parser')

	try: #try/catch in case a nucleotide/invalid sequence is entered
		url = soup.find(href=re.compile('http://raptorx.uchicago.edu/StructurePropertyPred/myjobs/')).get('href')

		zipid = url.split('_')
		zipid = zipid[1].replace('/','')

		raw = requests.get(url).text


		tree = html.fromstring(raw)

		treelist = tree.xpath('//*[@id="content"]/center[1]/text()')

		'''
		#No cancel
		while treelist != []:
			print('RaptorX Not Ready')
			time.sleep(20)
			raw = requests.get(url).text
			tree = html.fromstring(raw)
			treelist = tree.xpath('//*[@id="content"]/center[1]/text()')
		'''
		
		#Cancel after 20 min
		stime = time.time()
		while treelist != [] and time.time() < stime + 1200:
			print('RaptorX Not Ready')
			time.sleep(20)
			raw = requests.get(url).text
			tree = html.fromstring(raw)
			treelist = tree.xpath('//*[@id="content"]/center[1]/text()')

		if treelist == []:
			treelist = tree.xpath('//*[@id="infoval"]/script/text()')

			zippath = 'http://raptorx.uchicago.edu/'

			for word in treelist[0].split('"'):
				if word[:30] == '/StructurePropertyPred/myjobs/':
					zippath += word


			zipss = requests.get(zippath)

			while not zipss.ok:
				time.sleep(20)
				print('RaptorX Zip Not Ready')
			#zipss = batchtools.requestWait(zippath, 'RaptorX Zip Not Ready')
			
			mzip = request.urlopen(zippath)


			z= bytes(mzip.read())

			zf = ZipFile(io.BytesIO(z), "r")

			ss3path = zipid + '/' + zipid + '.ss3.txt' 

			sslist = zf.read(ss3path).decode('utf-8').splitlines()[2:]


			for i in range(len(sslist)):
				word = sslist[i].split()
				SS.pred += word[2]
				SS.hconf.append(word[3])
				SS.econf.append(word[4])
				SS.cconf.append(word[5])
				
			#SS.conf = "No conf given, Need formula to determine from h/e/c conf"	

			for i in range(len(SS.pred)):
				if SS.pred[i] == 'C':
					if SS.cconf[i] == "1.000":
						SS.conf += '9'			#puts in 9 if 100% confidence level
					else:
						SS.conf += SS.cconf[i][2]   #puts in the 3rd character from the individual conf value string e.g. 8 if string is 0.873
				if SS.pred[i] == 'E':
					if SS.econf[i] == "1.000":
						SS.conf += '9'
					else:
						SS.conf += SS.econf[i][2]
				if SS.pred[i] == 'H':
					if SS.hconf[i] == "1.000":
						SS.conf += '9'
					else:
						SS.conf += SS.hconf[i][2]

			SS.status = 1
			print("RaptorX Complete")
		else:
			SS.pred += "failed to respond in time"
			SS.conf += "failed to respond in time"
			SS.status = 2 #error status
			print("RaptorX failed: No response")
	except:
		SS.pred += "sequence not accepted"
		SS.conf += "sequence not accepted"
		SS.status = 4
		print("RaptorX failed: sequence not accepted")
		
	print("RAPTOR::")
	print(SS.pred)
	print(SS.conf)
	
	return SS