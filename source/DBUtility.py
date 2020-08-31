import subprocess
import pickle
import sqlite3
from sqlite3 import Error

def read_database_size(db_name, db_type):
    result = subprocess.run(["blastdbcmd", "-db", db_name, "-dbtype", db_type, "-info", "-exact_length"], stdout=subprocess.PIPE)
    db_info = result.stdout.decode('utf-8')
    db_info = db_info.split("\n")
    entries = db_info[1].split()
    num_seqs = int(entries[0].replace(',', ''))
    num_residues = int(entries[2].replace(',', ''))
    
    return num_seqs, num_residues

def load_db_parts(file_name):
    with open(file_name, "rb") as db_parts_file:
        db_parts = pickle.load(db_parts_file)
        
    return db_parts

def store_db_parts(db_parts, file_name):
    with open(file_name, "wb") as db_parts_file:
        pickle.dump(db_parts, db_parts_file)  


def read_database_parts(db_name, db_type):
    result = subprocess.run(["blastdbcmd", "-db", db_name, "-info", "-exact_length"], stdout=subprocess.PIPE)
    db_info = result.stdout.decode('utf-8').split("\n")
    relevant_lines = db_info[6:]
    num_lines = len(relevant_lines)
    parts = []
    for line in relevant_lines:
        part = line.strip().split("/")[-1]
        if part != "":
            parts.append(part)
    return parts

def create_database(db_file):
    """ create a database connection to a SQLite database """
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print(e)
    finally:
        conn.close()

# Create Database and Tables
def create_tables(database):
    query_string = """create table if not exists Query (
                                sequence text primary key
                                );"""

    db_instance_string = """create table if not exists DbInstance (
                                id integer primary key autoincrement, 
                                timestamp datetime default CURRENT_TIMESTAMP, 
                                dbname text, 
                                query text, 
                                partsfile text, 
                                foreign key(query) references Query(sequence)
                                );"""

    search_record_string = """create table if not exists SearchRecord (
                                dbinstanceid integer primary key,
                                record text,
                                foreign key(dbinstanceid) references DbInstance(id)
                                );"""
    try:
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        cursor.execute("PRAGMA foreign_keys = 1")
        cursor.execute(query_string)
        cursor.execute(db_instance_string)
        cursor.execute(search_record_string)
    except Error as e:
        print(e)
    finally:
        connection.close()


def initialize_database():
    create_database("incblast.db")
    create_tables("incblast.db")

def add_query(sequence):
    try:
        connection = sqlite3.connect("incblast.db")
        cursor = connection.cursor()
        cursor.execute("PRAGMA foreign_keys = 1")
        cursor.execute("insert into Query values(?)", (sequence,))
        connection.commit()
    except Error as e:
        print(e)
    finally:
        connection.close()

def add_db_instance(blast_db, query_sequence, parts_file):
    try:
        connection = sqlite3.connect("incblast.db")
        cursor = connection.cursor()
        cursor.execute("PRAGMA foreign_keys = 1")
        cursor.execute("insert into DbInstance(dbname, query, partsfile) values(?, ?, ?)",
                       (blast_db, query_sequence, parts_file))
        connection.commit()
    except Error as e:
        print(e)
    finally:
        connection.close()   

        
def add_search_record(db_instance_id, record_file):
    try:
        connection = sqlite3.connect("incblast.db")
        cursor = connection.cursor()
        cursor.execute("PRAGMA foreign_keys = 1")
        cursor.execute("insert into SearchRecord values(?, ?)", (db_instance_id, record_file))
        connection.commit()
    except Error as e:
        print(e)
    finally:
        connection.close()

def get_latest_dbinstance(query, database):
    try:
        connection = sqlite3.connect("incblast.db")
        cursor = connection.cursor()
        cursor.execute("""select id, timestamp, partsfile from DbInstance 
                            where id = (select max(id) from DbInstance where dbname = ? and query = ?)
                            """,
                       (database, query))
        result = cursor.fetchone()
        if result != None:
            dbinstanceid, timestamp, partsfile = result
            return dbinstanceid, timestamp, partsfile
    except Error as e:
        return None, None, None
    finally:
        connection.close() 
    return None, None, None

def execute_generic_command(command):
    command_tokens = command.split()
    result = subprocess.run(command_tokens, stdout=subprocess.PIPE)
    print(result.stdout.decode('utf-8'))
    return result.stdout.decode('utf-8')

def get_delta(blast_db, current_parts, db_type):
    whole_parts = read_database_parts(blast_db, db_type)
    whole_parts_set = set(whole_parts)
    current_parts_set = set(current_parts)
    delta_parts_set = whole_parts_set.difference(current_parts_set)

    if len(delta_parts_set) == 0:
    	return None, None
    
    print(current_parts_set)
    print(whole_parts_set)
    print(delta_parts_set)

    parts_string = ""
    for part in delta_parts_set:
        parts_string += part + " "
    parts_string = parts_string[:-1]
    
    delta_db = blast_db + "-delta"
    command = "blastdb_aliastool -dbtype " + db_type + " -dblist \""  + parts_string + "\" -out " + blast_db + "-delta -title " +  blast_db + "-delta"
    print("Command: ", command)
    execute_generic_command(command)
    
    return delta_db, delta_parts_set

def construct_command_string(commandOptions):
    program = commandOptions["program"]
    commandString = program
    for option in commandOptions:
        if option == "program":
            continue
        if option == "dbtype":
        	continue
        commandString += " -" + option + " " + str(commandOptions[option])
    print(commandString)
    return commandString

def parse_search_command(command):
    commandOptions = {}
    entries = command.split()
    commandOptions["program"] = entries[0]
    for i in range(1, len(entries), 2):
        option = entries[i][1:]
        value = entries[i+1]
        if value[0] == "-":
            print("Value was not provided for option \"" + option + "\"")
            return None
        commandOptions[option] = value.strip()

    print("program:", commandOptions["program"], ">")
    if commandOptions["program"] == "blastn":
    	commandOptions["dbtype"] = "nucl"
    	print("nucl search")
    else:
    	commandOptions["dbtype"] = "prot"
    	print("protein search")

    return commandOptions


def read_and_pickle_database_parts(db_name, db_type):
    result = subprocess.run(["blastdbcmd", "-db", db_name, "-dbtype", db_type, "-info", "-exact_length"], stdout=subprocess.PIPE)
    db_info = result.stdout.decode('utf-8').split("\n")
    relevant_lines = db_info[6:]
    num_lines = len(relevant_lines)
    parts = []
    for line in relevant_lines:
        part = line.strip().split("/")[-1]
        if part != "":
            parts.append(part)
            
    pickle_file_name = db_name + "-" + "parts.pkl"
    with open(pickle_file_name, "wb") as pickle_file:
        pickle.dump(parts, pickle_file)
        
    return parts, pickle_file_name

def get_search_record(dbinstanceid):
    try:
        connection = sqlite3.connect("incblast.db")
        cursor = connection.cursor()
        cursor.execute("""select record from SearchRecord 
                            where dbinstanceid = ?
                            """,
                       (dbinstanceid,))
        record = cursor.fetchone()
        if record != None:
            return record[0]
    except Error as e:
        return None
    finally:
        connection.close() 
    return None