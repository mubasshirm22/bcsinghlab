import os
import tempfile

from PSIPREDauto.functions import single_submit


def main():
    # Simple test sequence
    seq = (
        "PVRLIMAEHKRTWFKAVIIDNRSDQDHKKCGVEDTRRCWRHGRYGFMMEQHPISY"
        "CPANCYKDHANPMNEFQKIMFQWRSHMLVTTGGKKFMHTAAMLCFRELQMGKCN"
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "query.fasta")
        output_root = os.path.join(tmpdir, "psipred_out")

        # Write FASTA
        with open(fasta_path, "w") as f:
            f.write(">QUERY\n")
            f.write(seq + "\n")

        print("Submitting PSIPRED job via PSIPREDauto.single_submit...")
        # This will block until results are ready and write them under
        #   <output_root>/<fasta_file> output/
        single_submit(fasta_path, "no-reply@example.com", output_root, interval=1)

        # Locate .horiz file in the created output directory
        # single_submit creates a directory named "<fasta_file> output"
        out_dir_name = fasta_path + " output"
        horiz_path = None
        for root, dirs, files in os.walk(output_root):
            for fn in files:
                if fn.endswith(".horiz"):
                    horiz_path = os.path.join(root, fn)
                    break
            if horiz_path:
                break

        if not horiz_path or not os.path.exists(horiz_path):
            print("No .horiz file found in PSIPRED output.")
            return

        print("\nPreview of .horiz file:", horiz_path)
        with open(horiz_path, "r") as hf:
            for line in hf:
                line = line.strip()
                if line.startswith("Conf") or line.startswith("Pred"):
                    print(line)


if __name__ == "__main__":
    main()