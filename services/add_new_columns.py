import os
import psycopg2
from psycopg2 import sql

DATABASE_URL = os.environ.get('DATABASE_URL')

if not DATABASE_URL:
	print("ERROR: DATABASE_URL environment variable not set")
	exit(1)

conn = psycopg2.connect(DATABASE_URL)
cursor = conn.cursor()

# Add new columns for PHDpsi, PROFsec, and Predator if they don't exist
new_columns = [
	('PHDPSIPRED', 'TEXT'),
	('PHDPSICONF', 'TEXT'),
	('PHDPSISTAT', 'INT'),
	('PHDPSIMSG', 'TEXT'),
	('PROFSECPRED', 'TEXT'),
	('PROFSECCONF', 'TEXT'),
	('PROFSECSTAT', 'INT'),
	('PROFSECMSG', 'TEXT'),
	('PREDATORPRED', 'TEXT'),
	('PREDATORCONF', 'TEXT'),
	('PREDATORSTAT', 'INT'),
	('PREDATORMSG', 'TEXT'),
]

for col_name, col_type in new_columns:
	try:
		# Check if column exists
		cursor.execute("""
			SELECT column_name 
			FROM information_schema.columns 
			WHERE table_name='seqtable' AND column_name=%s
		""", (col_name.lower(),))
		
		if cursor.fetchone() is None:
			# Column doesn't exist, add it
			cursor.execute(
				sql.SQL("ALTER TABLE seqtable ADD COLUMN {} {}")
				.format(sql.Identifier(col_name.lower()), sql.SQL(col_type))
			)
			print(f"Added column: {col_name}")
		else:
			print(f"Column {col_name} already exists")
	except Exception as e:
		print(f"Error adding column {col_name}: {e}")

conn.commit()
cursor.close()
conn.close()

print("Done! New columns added to database.")
