import DBUtility as DB
import iBlast

def create_first_blastdb(prefix, num_parts, db_type, output):
	command = "blastdb_aliastool -dblist "
	parts = "\""
	for i in range(num_parts):
		part = prefix + ".{:02d} ".format(i)
		parts += part
	parts = parts.strip()


	parts += "\" "
	command += parts
	command += "-dbtype " + db_type + " "
	command += "-out " + output + " "
	command += "-title " + output

	print(command)
	return command

DB.init_db()
#command = create_first_blastdb("../blastn/nt", 2, "nucl", "../blastn/ntsq") 



#iBlast.incremental_blast("blastn -db ../blastn/ntsq -query ../blastn/small-query.fasta -outfmt 5 -out ../blastn/sq-result.xml")
