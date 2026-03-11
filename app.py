import json
import time
import os
from lxml import html
from services import ss, psi, jpred, raptorx, pss, sable, sspro, yaspin, phdpsi, profsec, predator, emailtools, htmlmaker, batchtools, maketable
from datetime import datetime

from forms import SubmissionForm
from flask import Flask, render_template, request, current_app, send_file, redirect, url_for, jsonify
import threading
import secrets
import requests

import psycopg2
from psycopg2.extras import RealDictCursor
from psycopg2 import sql

DATABASE_URL = os.environ.get('DATABASE_URL')



def dbselect(rowid):
	conn = psycopg2.connect(DATABASE_URL)
	cursor = conn.cursor(cursor_factory=RealDictCursor)
	cursor.execute("SELECT * FROM seqtable WHERE ID = (%s)",(rowid,))
	jsonresults = json.dumps(cursor.fetchall(), indent=2)
	cursor.close()
	return jsonresults

def dbdelete():
	conn = psycopg2.connect(DATABASE_URL)
	cursor = conn.cursor()
	cursor.execute("SELECT COUNT(*) FROM seqtable")
	numrowsdb = cursor.fetchall()
	numrowsdb = numrowsdb[0][0]

	if numrowsdb > 8000:
		cursor.execute(
				'''
				DELETE FROM seqtable
				WHERE ID = any (array(SELECT ID FROM seqtable ORDER BY convert_to(ID, 'SQL_ASCII') ASC LIMIT 1000))
				''')
		conn.commit()
	cursor.close()

def dbinsert(rowid, rowseq):
	conn = psycopg2.connect(DATABASE_URL)
	dbdelete() #Deletes 1000 oldest rows if table is larger than 8000 rows
	cursor = conn.cursor()
	cursor.execute("INSERT INTO seqtable (ID, SEQ) VALUES (%s, %s)", (rowid, rowseq))
	conn.commit()
	cursor.close()

def dbupdate(rowid, rowcol, rowval):
	conn = psycopg2.connect(DATABASE_URL)
	cursor = conn.cursor()
	cursor.execute(
		sql.SQL("UPDATE seqtable SET {} = (%s) WHERE ID = (%s)")
		.format(sql.Identifier(rowcol.lower())),(rowval, rowid))
	conn.commit()
	cursor.close()

def check_service_health():
	"""Quick health check for all services"""
	service_urls = {
		"JPred": "http://www.compbio.dundee.ac.uk/jpred4/",
		"PSI": "http://bioinf.cs.ucl.ac.uk/psipred/",
		"Sable": "http://sable.cchmc.org/",
		"SSPro": "http://scratch.proteomics.ics.uci.edu/",
		"Yaspin": "http://www.ibi.vu.nl/programs/yaspinwww/",
		"PHDpsi": "https://predictprotein.org/",
		"PROFsec": "https://predictprotein.org/",
		"Predator": "http://npsa-pbil.ibcp.fr/"
	}

	status = {}
	for service_name, url in service_urls.items():
		try:
			response = requests.get(url, timeout=5)
			if response.status_code == 200:
				status[service_name] = "UP"
			else:
				status[service_name] = "DOWN"
		except:
			status[service_name] = "DOWN"

	# Check STRIDE (PDB service) - uses POST endpoint
	try:
		stride_url = "https://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py"
		# Try a simple POST with minimal data to check if server is reachable
		test_data = {
			"pdbid": "TEST",
			"paste_field": "",
			"action": "compute",
			"contact_threshold": "6",
			"sensitive": "true"
		}
		response = requests.post(stride_url, data=test_data, timeout=5)
		# Server should respond (even if with an error about invalid PDB ID)
		if response.status_code == 200:
			status["STRIDE"] = "UP"
		else:
			status["STRIDE"] = "DOWN"
	except:
		status["STRIDE"] = "DOWN"

	return status

#Dictionary containing sites and their classes
siteDict = {
	"JPred": jpred,
	"PSI": psi,
	#"PSS": pss,
	#"RaptorX": raptorx,
	"Sable": sable,
	"Yaspin": yaspin,
	"SSPro": sspro,
	"PHDpsi": phdpsi,
	"PROFsec": profsec,
	"Predator": predator
}

siteLimit = {
	"JPred": 20,
	"PSI": 20,
	#"PSS": 3,
	#"RaptorX": 20,
	"Sable": 20,
	"Yaspin": 3,
	"SSPro": 5,
	"PHDpsi": 5,
	"PROFsec": 5,
	"Predator": 5
}


app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET')
if app.config['SECRET_KEY'] is None:
	app.config['SECRET_KEY'] = secrets.token_urlsafe(16)

#Login to email account to be able to send emails
#email_service = emailtools.login()
#email = emailtools.getEmailAddress(email_service)

#Url of hosted site
siteurl = os.environ.get('SITE_URL')
if siteurl is None :
	siteurl = ""


# ---------------------------------------------------------------------------
# Lab site routes
# ---------------------------------------------------------------------------

@app.route('/')
def lab_home():
	return render_template('lab/home.html')

@app.route('/research')
def lab_research():
	return render_template('lab/research.html')

@app.route('/team')
def lab_team():
	return render_template('lab/team.html')

@app.route('/tools')
def lab_tools():
	return render_template('lab/tools.html')

@app.route('/tutorials')
def lab_tutorials():
	return render_template('lab/tutorials.html')

@app.route('/contact')
def lab_contact():
	return render_template('lab/contact.html')


# ---------------------------------------------------------------------------
# SSPred tool routes (moved from / to /tools/sspred)
# ---------------------------------------------------------------------------

@app.route('/tools/sspred', methods = ['GET', 'POST'])
def sspred_index(name=None):
	form = SubmissionForm()
	print(threading.activeCount())
	runningCounter = {
		"JPred": 0,
		"PSI": 0,
		#"PSS": 0,
		#"RaptorX": 0,
		"Sable": 0,
		"Yaspin": 0,
		"SSPro": 0,
		"PHDpsi": 0,
		"PROFsec": 0,
		"Predator": 0
	}
	for t in threading.enumerate():
		if t.getName() in runningCounter.keys():
			runningCounter[t.getName()] += 1

	# Check service health
	serviceStatus = check_service_health()

	if form.validate_on_submit():

		if threading.activeCount() > 100:
			return redirect(url_for('sspred_error'))
		# Normalize the typed sequence (remove whitespace) up front
		clean_seq = ''.join(form.seqtext.data.split())
		print("DEBUG: original form sequence length:", len(form.seqtext.data or ""))
		print("DEBUG: cleaned form sequence length:", len(clean_seq))

		post_data = {
			'seqtext': clean_seq,
			'email': form.email.data,
			'JPred': form.JPred.data,
			'PSI':   form.PSI.data,
			#'PSS':   form.PSS.data,
			#'RaptorX': form.RaptorX.data,
			'Sable':   form.Sable.data,
			'Yaspin':   form.Yaspin.data,
			'SSPro':   form.SSPro.data,
			'PHDpsi':   form.PHDpsi.data,
			'PROFsec':   form.PROFsec.data,
			'Predator':   form.Predator.data,
			'submitbtn': 'Submit'
			}

		total_sites = validate_sites(post_data)
		post_data.update({'total_sites' : total_sites, 'completed': 0}) # add total into to post_data dictionary and a completed prediction counter
		print(post_data)
		seq = post_data['seqtext']
		print("DEBUG: initial seq length used for prediction:", len(seq))

		startTime = batchtools.randBase62()

		#Stores currently completed predictions
		ssObject = []

		# Insert initial sequence (may be empty; will be updated if PDB succeeds)
		dbinsert(startTime, seq)

		pdbdata = None
		if form.structureId.data is not None:
			# Auto-convert PDB ID to uppercase (PDB IDs are always uppercase)
			pdb_id = form.structureId.data.upper().strip()
			chain_id = form.chainId.data.upper().strip() if form.chainId.data else None
			pdbdata = batchtools.pdbget(pdb_id, chain_id)
			if pdbdata is not None:
				dbupdate(startTime, 'pdb', json.dumps(pdbdata))
				dbupdate(startTime, 'seq', pdbdata['primary'])
				seq = pdbdata['primary']
				print("DEBUG: using PDB-derived sequence, length:", len(seq))
			else:
				print("DEBUG: pdbget returned None for", pdb_id, chain_id)

		# If after PDB lookup we still have no sequence, abort cleanly
		if not seq:
			statusMsg = "No valid sequence found. Please provide a sequence or a valid structure ID and chain."
			return redirect(url_for('sspred_dboutput', var=startTime))

		sendData(seq, startTime, ssObject, post_data, pdbdata)
		return redirect(url_for('sspred_dboutput', var=startTime))

	return render_template('index.html', form=form, counter=runningCounter, serviceStatus=serviceStatus)

@app.route('/tools/sspred/error/')
def sspred_error():
	return('There are too many jobs running, please try again later')

@app.route('/tools/sspred/archive/<page>')
def sspred_archive(page):
	if page[0] == '0':
		return("Page not found")
	if page.isdigit():
		if int(page) >= 1:
			namelist = []
			timelist= []
			seqlist= []
			conn = psycopg2.connect(DATABASE_URL)
			cursor = conn.cursor(cursor_factory=RealDictCursor)
			limit = 20
			offset = int(page) -1
			offset = offset * limit
			cursor.execute('''
					SELECT id, seq
					FROM seqtable
					ORDER BY ID DESC LIMIT %s OFFSET %s
			''',(limit, offset))
			jsonresults = json.dumps(cursor.fetchall(), indent=2)

			cursor.close()

			return render_template('archives.html', data=jsonresults, pagenum=page)
	else:
		return("Page not found")

@app.route('/tools/sspred/archive')
def sspred_archive_redirect():
	return redirect(url_for('sspred_archive', page=1))

@app.route('/tools/sspred/output/<var>')
def sspred_output(var):
	print('output/'+var+'/'+var+'.html')
	try:
		return send_file('output/'+var+'/'+var+'.html')
	except Exception as e:
		return "not found"

@app.route('/tools/sspred/dboutput/<var>')
def sspred_dboutput(var):
	outputjson = dbselect(var)
	if outputjson == "[]":
		return "not found"
	try:
		if request.args.get('json') == '1':
			return outputjson
		return render_template('dboutput.html', data=outputjson)
	except Exception as e:
		return "not found"


# ---------------------------------------------------------------------------
# Legacy 301 redirects (backward compatibility for bookmarked/shared links)
# ---------------------------------------------------------------------------

@app.route('/archive')
def legacy_archive():
	return redirect(url_for('sspred_archive_redirect'), 301)

@app.route('/archive/<page>')
def legacy_showall(page):
	return redirect(url_for('sspred_archive', page=page), 301)

@app.route('/output/<var>')
def legacy_output(var):
	return redirect(url_for('sspred_output', var=var), 301)

@app.route('/dboutput/<var>')
def legacy_dboutput(var):
	return redirect(url_for('sspred_dboutput', var=var), 301)

@app.route('/error/')
def legacy_error():
	return redirect(url_for('sspred_error'), 301)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def run(predService, seq, name, ssObject,
 startTime, post_data, pdbdata):
	try:
		tcount = 0
		for t in threading.enumerate():
			if t.getName() == name:
				tcount += 1

		if tcount > siteLimit[name]:
			tempSS = ss.SS(name)
			tempSS.pred = "Queue Full"
			tempSS.conf = "Queue Full"
			tempSS.status = -1
			dbupdate(startTime, name + "msg", name + " – didn't run because the remote queue is full.")
		else:
			import time as time_module
			start_time = time_module.time()
			dbupdate(startTime, name + "msg", name + " – submitting job to remote server...")
			try:
				#tempSS = predService.get(seq, tcount)
				tempSS = predService.get(seq)
				elapsed_min = int((time_module.time() - start_time) / 60)
				if tempSS.status == 1 or tempSS.status == 3:
					dbupdate(startTime, name + "msg",
					         name + " – finished successfully in " + str(elapsed_min) + " minutes.")
				elif tempSS.status == 2:
					if "failed to respond after" in tempSS.pred:
						dbupdate(startTime, name + "msg",
						         name + " – gave up after " + str(elapsed_min) +
						         " minutes without a response (remote server might be slow or down).")
					elif "Queue Full" in tempSS.pred:
						dbupdate(startTime, name + "msg",
						         name + " – didn't run because the remote queue is full.")
					else:
						dbupdate(startTime, name + "msg",
						         name + " – failed: " + tempSS.pred)
				elif tempSS.status == 4:
					dbupdate(startTime, name + "msg",
					         name + " – sequence not accepted by the server.")
				else:
					dbupdate(startTime, name + "msg",
					         name + " – still waiting after " + str(elapsed_min) +
					         " minutes (may give up if the server doesn't respond).")
			except Exception as e:
				tempSS = ss.SS(name)
				tempSS.pred = "Service Error: " + str(e)
				tempSS.conf = "Service Error: " + str(e)
				tempSS.status = 2
				dbupdate(startTime, name + "msg", name + " – didn't run because of a service error: " + str(e))
				print(name + " failed with exception: " + str(e))

		dbupdate(startTime, tempSS.name + "pred", tempSS.pred)
		dbupdate(startTime, tempSS.name + "conf", tempSS.conf)
		dbupdate(startTime, tempSS.name + "stat", tempSS.status)

		ssObject.append(tempSS)
		majority = batchtools.majorityVote(seq, ssObject)
		dbupdate(startTime, 'majorityvote', majority)

		post_data['completed'] += 1
		if post_data['completed'] == post_data['total_sites']:
			print("All predictions completed.")
			statusMsg = "All services complete."
			failedServices = []
			for ssobj in ssObject:
				if ssobj.status != 1 and ssobj.status != 3:
					if ssobj.status == -1:
						failedServices.append(ssobj.name + " didn't run because queue is full")
					elif ssobj.status == 2:
						if "Service Error" in ssobj.pred:
							failedServices.append(ssobj.name + " didn't run because " + ssobj.pred.replace("Service Error: ", ""))
						else:
							failedServices.append(ssobj.name + " didn't run because " + ssobj.pred)
					elif ssobj.status == 4:
						failedServices.append(ssobj.name + " didn't run because sequence not accepted")
					else:
						failedServices.append(ssobj.name + " didn't complete")
			if failedServices:
				statusMsg += " " + ". ".join(failedServices) + "."
			# dbupdate(startTime, 'status', statusMsg)  # Commented out: status column doesn't exist in database
			if post_data['email'] != "": #if all completed and user email is not empty, send email
				print("Sending results notification to " + post_data['email'])
				emailtools.send_job_notification(post_data['email'], startTime, siteurl)
	except Exception as e:
		# Catch any unexpected errors to ensure thread completes
		print(name + " thread failed with unexpected error: " + str(e))
		import traceback
		traceback.print_exc()
		try:
			tempSS = ss.SS(name)
			tempSS.pred = "Thread Error: " + str(e)
			tempSS.conf = "Thread Error: " + str(e)
			tempSS.status = 2
			dbupdate(startTime, name + "msg", name + " didn't run because thread error: " + str(e))
			dbupdate(startTime, tempSS.name + "pred", tempSS.pred)
			dbupdate(startTime, tempSS.name + "conf", tempSS.conf)
			dbupdate(startTime, tempSS.name + "stat", tempSS.status)
			# Check if this service is already in ssObject by name
			found = False
			for ssobj in ssObject:
				if ssobj.name == name:
					found = True
					break
			if not found:
				ssObject.append(tempSS)
			post_data['completed'] += 1
			if post_data['completed'] == post_data['total_sites']:
				statusMsg = "All services complete. " + name + " didn't run because thread error: " + str(e) + "."
				# dbupdate(startTime, 'status', statusMsg)  # Commented out: status column doesn't exist in database
		except Exception as e2:
			print("Failed to update database after thread error for " + name + ": " + str(e2))
			import traceback
			traceback.print_exc()


#Sends sequence based off whatever was selected before submission
def sendData(seq, startTime, ssObject, post_data, pdbdata):
	for key in post_data.keys():
		if key in siteDict:
			if post_data[key]:
				mythread = threading.Thread(target=run, args=(siteDict[key], seq, key, ssObject, startTime, post_data, pdbdata))
				mythread.setName(key)
				mythread.start()
				print("Sending sequence to " + key)

#Takes a form from post and returns the number of sites selected.
def validate_sites(form):
	count = 0
	for key in siteDict.keys():
		if form[key]:
			count += 1
	return count

# ---------------------------------------------------------------------------
# ProtPipe — protein analysis pipeline routes
# ---------------------------------------------------------------------------

try:
	from pipeline.runner import submit_job as _pipe_submit, get_status as _pipe_status, get_summary as _pipe_summary
	_PIPELINE_AVAILABLE = True
except ImportError as _pipe_err:
	_PIPELINE_AVAILABLE = False
	print(f"[protpipe] pipeline not available: {_pipe_err}")


@app.route('/tools/protpipe', methods=['GET', 'POST'])
def protpipe_index():
	if not _PIPELINE_AVAILABLE:
		return render_template('protpipe/unavailable.html',
			reason="Pipeline dependencies not installed. Run: pip install biopython svgwrite"), 503

	if request.method == 'POST':
		input_data = {
			'input_type':      request.form.get('input_type', 'raw_fasta'),
			'sequence_input':  request.form.get('sequence_input', '').strip(),
			'run_blast':       'run_blast'   in request.form,
			'run_hmmer':       'run_hmmer'   in request.form,
			'run_phobius':     'run_phobius' in request.form,
		}
		if not input_data['sequence_input']:
			return render_template('protpipe/index.html', error="Please enter a sequence or accession.")

		job_id = _pipe_submit(input_data)
		return redirect(url_for('protpipe_results', job_id=job_id))

	return render_template('protpipe/index.html', error=None)


@app.route('/tools/protpipe/results/<job_id>')
def protpipe_results(job_id):
	if not _PIPELINE_AVAILABLE:
		return render_template('protpipe/unavailable.html',
			reason="Pipeline dependencies not installed."), 503

	status_data = _pipe_status(job_id)
	if status_data.get('status') == 'not_found':
		return render_template('protpipe/index.html', error=f"Job '{job_id}' not found."), 404

	summary = {}
	if status_data.get('status') == 'complete':
		summary = _pipe_summary(job_id)

	return render_template('protpipe/results.html',
		job_id=job_id,
		status_data=status_data,
		summary=summary,
	)


@app.route('/tools/protpipe/status/<job_id>')
def protpipe_status_api(job_id):
	"""JSON polling endpoint — called by the results page every few seconds."""
	return jsonify(_pipe_status(job_id) if _PIPELINE_AVAILABLE else {"status": "error"})


@app.route('/tools/protpipe/figure/<job_id>')
def protpipe_figure(job_id):
	"""Serve the SVG domain figure for a completed job."""
	import os
	from pipeline.utils.jobs import job_dir
	svg_path = os.path.join(job_dir(job_id), "domain_figure.svg")
	if not os.path.exists(svg_path):
		return ("Figure not available", 404)
	return send_file(svg_path, mimetype="image/svg+xml")


@app.route('/tools/protpipe/archive')
def protpipe_archive():
	return render_template('protpipe/archive.html')


@app.route('/tools/protpipe/download/<job_id>/json')
def protpipe_download_json(job_id):
	if not _PIPELINE_AVAILABLE:
		return ("Pipeline not available", 503)
	from pipeline.utils.jobs import job_dir, get_summary
	import io
	summary = get_summary(job_id)
	if not summary:
		return ("Job not found", 404)
	data = json.dumps(summary, indent=2)
	return send_file(
		io.BytesIO(data.encode()),
		mimetype="application/json",
		as_attachment=True,
		attachment_filename=f"protpipe_{job_id}.json",
	)


@app.route('/tools/protpipe/download/<job_id>/svg')
def protpipe_download_svg(job_id):
	if not _PIPELINE_AVAILABLE:
		return ("Pipeline not available", 503)
	from pipeline.utils.jobs import job_dir as _job_dir
	svg_path = os.path.join(_job_dir(job_id), "domain_figure.svg")
	if not os.path.exists(svg_path):
		return ("Figure not available", 404)
	return send_file(
		svg_path,
		mimetype="image/svg+xml",
		as_attachment=True,
		attachment_filename=f"protpipe_{job_id}_figure.svg",
	)


@app.route('/tools/protpipe/download/<job_id>/txt')
def protpipe_download_txt(job_id):
	if not _PIPELINE_AVAILABLE:
		return ("Pipeline not available", 503)
	from pipeline.utils.jobs import get_summary
	import io
	summary = get_summary(job_id)
	if not summary:
		return ("Job not found", 404)
	lines = _make_text_report(job_id, summary)
	return send_file(
		io.BytesIO(lines.encode()),
		mimetype="text/plain",
		as_attachment=True,
		attachment_filename=f"protpipe_{job_id}.txt",
	)


def _make_text_report(job_id, summary):
	"""Generate a plain-text summary report from a completed pipeline job.

	Note: summary["properties"], summary["blast"], summary["hmmer"], and
	summary["phobius"] are the raw data dicts (no outer "status" wrapper) —
	that unwrapping happens in runner.py when summary.json is assembled.
	summary["retrieval"] IS the full module result including "status".
	"""
	lines = []
	lines.append("=" * 60)
	lines.append("ProtPipe — Protein Analysis Report")
	lines.append(f"Job ID : {job_id}")
	lines.append(f"Date   : {summary.get('completed_at', 'N/A')}")
	lines.append("=" * 60)

	# Retrieval (full result dict with status key)
	ret = summary.get("retrieval", {})
	if ret.get("status") == "ok":
		lines.append(f"\nSequence Source : {ret.get('source','N/A')}")
		lines.append(f"Header          : {ret.get('header','N/A')}")
		lines.append(f"Organism        : {ret.get('organism','N/A')}")
		seq = ret.get("sequence","")
		lines.append(f"Length          : {len(seq)} aa")
		if seq:
			lines.append(f"\nSequence:\n{seq}")

	# Properties (data dict directly — no status key)
	prop = summary.get("properties", {})
	if prop:
		lines.append("\n--- Physicochemical Properties ---")
		if prop.get('molecular_weight_da'): lines.append(f"Molecular Weight : {prop['molecular_weight_da']} Da")
		if prop.get('isoelectric_point'):   lines.append(f"Isoelectric Point: {prop['isoelectric_point']}")
		if prop.get('gravy') is not None:   lines.append(f"GRAVY            : {prop['gravy']}")
		if prop.get('instability_index') is not None: lines.append(f"Instability Index: {prop['instability_index']}")
		if prop.get('aromaticity') is not None: lines.append(f"Aromaticity      : {prop['aromaticity']}")

	# BLAST (data dict with "hits" key directly)
	hits = summary.get("blast", {}).get("hits", [])
	if hits:
		lines.append(f"\n--- BLAST Homology ({len(hits)} hits) ---")
		for i, h in enumerate(hits, 1):
			lines.append(
				f"{i:2}. {h.get('accession','?')}  {h.get('title','?')[:50]}"
				f"  E={h.get('e_value','?')}  Id={h.get('identity_pct','?')}%  Cov={h.get('coverage_pct','?')}%"
			)

	# HMMER (data dict with "domains" key directly)
	domains = summary.get("hmmer", {}).get("domains", [])
	if domains:
		lines.append(f"\n--- Pfam Domains ({len(domains)} found) ---")
		for d in domains:
			lines.append(
				f"  {d.get('name','?')}  [{d.get('seq_start','?')}-{d.get('seq_end','?')}]"
				f"  E={d.get('e_value','?')}  {d.get('description','')}"
			)

	# Phobius (data dict with topology keys directly)
	phob = summary.get("phobius", {})
	if phob:
		lines.append("\n--- Signal Peptide & TM Topology ---")
		sp = "Yes, cleavage after position " + str(phob.get('signal_peptide_end')) if phob.get('has_signal_peptide') else "Not detected"
		lines.append(f"Signal Peptide : {sp}")
		lines.append(f"TM Helices     : {phob.get('tm_count', 0)}")
		if phob.get('topology'):
			lines.append(f"Topology       : {phob['topology']}")

	lines.append("\n" + "=" * 60)
	lines.append("Generated by ProtPipe — Singh Lab, Brooklyn College")
	lines.append("=" * 60)
	return "\n".join(lines)


if __name__ == "__main__":
	app.run(debug=True) #Run on localhost 127.0.0.1:5000
	#app.run(host='0.0.0.0', debug=True) #Run online on public IP:5000
