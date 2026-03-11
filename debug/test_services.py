import os
import sys
import time

import requests
from guerrillamail import GuerrillaMailSession

"""
Small standalone script to debug external secondary-structure services
without going through the Flask app.

Usage (from project root):

    source venv/bin/activate
    python debug/test_services.py predator
    python debug/test_services.py yaspin
    python debug/test_services.py psi
    python debug/test_services.py phdpsi
    python debug/test_services.py profsec

Edit TEST_SEQ below to use your own sequence.
This script just prints raw HTTP status, headers and the first chunk
of response text to your terminal so we can see exactly what each
remote server is doing.
"""


# *** CHANGE THIS SEQUENCE WHEN TESTING ***
TEST_SEQ = (
    "MKTAYIAKQRQISFVKSHFSRQDILDLICRHLQQNTPANVVVTDDDLRDRLA"
    "KMAGLGTDEEVISVSEYQGQ"
)


def _print_response(label: str, r: requests.Response, max_chars: int | None = None) -> None:
    """Utility: pretty‑print HTTP response info.

    If max_chars is None, print the full body. Otherwise truncate.
    """
    print(f"\n==== {label} RESPONSE ====")
    print("URL:", r.url)
    print("Status code:", r.status_code)
    print("Headers:")
    for k, v in r.headers.items():
        print(f"  {k}: {v}")
    text = r.text or ""
    if max_chars is None:
        print("\nBody (full):")
        print(text)
    else:
        print(f"\nBody (first {max_chars} chars):")
        print(text[:max_chars])
    print("==== END RESPONSE ====\n")


def test_predator(seq: str) -> None:
    """
    Debug PREDATOR by calling the same endpoint the wrapper uses.
    This prints the FULL HTML of whatever page it returns so we can
    scroll around and see exactly where (or if) the prediction appears.
    """
    url = "https://npsa-prabi.ibcp.fr/cgi-bin/secpred_preda.pl"
    payload = {
        "title": "",
        "notice": seq,
        "ali_width": "70",
        "predatorssmat": "dssp",
    }
    print("Testing Predator...")
    print("POST", url)
    try:
        r = requests.post(
            url,
            data=payload,
            timeout=180,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
    except Exception as e:
        print("Predator request failed with exception:", repr(e))
        return
    # Print the entire HTML body so you can scroll and search for H/E/C lines
    _print_response("Predator", r, max_chars=None)


def _poll_url(url: str, label: str, interval: int = 20, timeout: int = 2700) -> None:
    """
    Simple polling helper: repeatedly GET a URL and dump the first successful
    response (or print a timeout message).
    """
    print(f"\nPolling {label} URL:", url)
    start = time.time()
    attempt = 0
    while time.time() - start < timeout:
        attempt += 1
        try:
            r = requests.get(url, timeout=30)
        except Exception as e:
            print(f"[{label}] attempt {attempt}: request error:", repr(e))
            time.sleep(interval)
            continue

        # For debugging we just want the first non‑empty response;
        # some services return an HTML shell first and then full results.
        if r.ok and r.text and r.text.strip():
            print(f"[{label}] got non‑empty response on attempt {attempt}")
            _print_response(label, r, max_chars=None)
            return

        print(f"[{label}] attempt {attempt}: still not ready (status {r.status_code})")
        time.sleep(interval)

    print(f"[{label}] gave up after {int((time.time() - start)/60)} minutes without usable response.")


def test_yaspin(seq: str) -> None:
    """
    Debug YASPIN: submit a job and then poll the results.out URL,
    dumping whatever comes back.
    """
    print("Testing Yaspin...")
    session = GuerrillaMailSession()
    email_address = session.get_session_state()["email_address"]

    payload = {
        "seq": seq,
        "mbjob[description]": "testprot",
        "nnmethod": "dssp",
        "smethod": "nr",
        "yaspin_align": "YASPIN prediction",
        "email": email_address,
    }
    files = {"seq_file": ""}

    submit_url = "http://www.ibi.vu.nl/programs/yaspinwww/"
    print("POST", submit_url)
    try:
        r = requests.post(submit_url, data=payload, files=files, timeout=60)
    except Exception as e:
        print("Yaspin submit failed with exception:", repr(e))
        return

    _print_response("Yaspin submit", r)
    result_url = r.url + "results.out"
    _poll_url(result_url, "Yaspin results.out")


def test_psipred(seq: str) -> None:
    """
    Debug PSIPRED (PSI) by reproducing the API calls from psi.py.
    Note: at the moment their SSL certificate is expired, so this
    will likely fail with a certificate error – which is what we want to see.
    """
    print("Testing PSIPRED (PSI)...")
    session = GuerrillaMailSession()
    email_address = session.get_session_state()["email_address"]

    url = "http://bioinf.cs.ucl.ac.uk/psipred/api/submission/"
    fasta_seq = ">tr|B0R5N9|B0R5N9\n" + seq
    payload = {"input_data": fasta_seq}
    data = {
        "job": "psipred",
        "submission_name": "testing",
        "email": email_address,
    }

    print("POST", url)
    try:
        r = requests.post(url, data=data, files=payload, headers={"accept": "application/json"}, timeout=30)
    except Exception as e:
        print("PSIPRED submit failed with exception:", repr(e))
        return
    _print_response("PSIPRED submit", r)


def test_phdpsi(seq: str) -> None:
    """
    Debug PHDpsi via the approximate PredictProtein endpoint.
    This just prints whatever HTML/JSON they send back.
    """
    print("Testing PHDpsi (PredictProtein)...")
    session = GuerrillaMailSession()
    email_address = session.get_session_state()["email_address"]

    url = "https://predictprotein.org/api/submit"
    payload = {
        "sequence": seq,
        "email": email_address,
        "jobname": "phdpsi_test",
        "output": "html",
    }

    print("POST", url)
    try:
        r = requests.post(url, data=payload, timeout=60)
    except Exception as e:
        print("PHDpsi submit failed with exception:", repr(e))
        return
    _print_response("PHDpsi submit", r)


def test_profsec(seq: str) -> None:
    """
    Debug PROFsec via the same PredictProtein endpoint.
    """
    print("Testing PROFsec (PredictProtein)...")
    session = GuerrillaMailSession()
    email_address = session.get_session_state()["email_address"]

    url = "https://predictprotein.org/api/submit"
    payload = {
        "sequence": seq,
        "email": email_address,
        "jobname": "profsec_test",
        "output": "html",
    }

    print("POST", url)
    try:
        r = requests.post(url, data=payload, timeout=60)
    except Exception as e:
        print("PROFsec submit failed with exception:", repr(e))
        return
    _print_response("PROFsec submit", r)


def main():
    if len(sys.argv) < 2:
        print(
            "Usage: python debug/test_services.py "
            "[predator|yaspin|psi|phdpsi|profsec]\n"
            "Edit TEST_SEQ in this file to change the test sequence."
        )
        return

    service = sys.argv[1].lower()
    seq = TEST_SEQ

    if service == "predator":
        test_predator(seq)
    elif service == "yaspin":
        test_yaspin(seq)
    elif service == "psi" or service == "psipred":
        test_psipred(seq)
    elif service == "phdpsi":
        test_phdpsi(seq)
    elif service == "profsec":
        test_profsec(seq)
    else:
        print("Unknown service:", service)
        print("Valid options: predator, yaspin, psi, psipred, phdpsi, profsec")


if __name__ == "__main__":
    main()

