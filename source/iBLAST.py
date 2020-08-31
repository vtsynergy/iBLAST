import pickle
import sqlite3
from collections import namedtuple
import DBUtility as DB
import subprocess
from sqlite3 import Error
from BlastnMergerModule import BlastnMerger
from BlastpMergerModule import BlastpMerger
import sys
import copy

def perform_fresh_search(command):
    command_options = DB.parse_search_command(command)
    blast_db = command_options["db"]
    query_sequence = command_options["query"]
    record_file = command_options["out"]
    db_parts, parts_file = DB.read_and_pickle_database_parts(command_options["db"], command_options["dbtype"])

    
    output = DB.execute_generic_command(command)
    print(output)
    DB.add_query(query_sequence)
    DB.add_db_instance(blast_db, query_sequence, parts_file)
    dbinstanceid, timestamp, partsfile = DB.get_latest_dbinstance(query_sequence, blast_db)
    print(dbinstanceid)
    DB.add_search_record(dbinstanceid, record_file)


def perform_incrermental_search(command, command_options, dbinstanceid, current_partsfile):
    print("output file: ", command_options["out"])
    blast_db = command_options["db"]
    query_sequence = command_options["query"]
    record_file = command_options["out"]
    current_parts = DB.load_db_parts(current_partsfile)    
    delta_db, delta_parts = DB.get_delta(blast_db, current_parts, command_options["dbtype"])
    if(delta_db == None):
    	print("The database has not changed since last search. Displaying result from last search.")
    	return
    
    delta_command_options = copy.deepcopy(command_options)
    delta_command_options["db"] = delta_db
    delta_command_options["out"] = command_options["out"] + "-delta.xml"
    delta_command = DB.construct_command_string(delta_command_options)
    DB.execute_generic_command(delta_command)
    print("output file: ", command_options["out"])
    current_result = DB.get_search_record(dbinstanceid)
    delta_result = delta_command_options["out"]
    if command_options["dbtype"] == "nucl":
        merger = BlastnMerger()
        #current, delta, merged = 
        merger.merge_multi(current_result, delta_result, command_options["out"])
    if command_options["dbtype"] == "prot":
        merger = BlastpMerger()
        #current, delta, merged = 
        merger.merge_multi(current_result, delta_result, command_options["out"])


def incremental_blast(command):
    commandOptions = DB.parse_search_command(command)
    dbname = commandOptions["db"]
    query = commandOptions["query"]
    print(dbname, query)
    dbinstanceid, timestamp, current_partsfile = DB.get_latest_dbinstance(query, dbname)
    
    if(dbinstanceid == None):
        print("No record found for this query. A fresh blast search will be performed")
        perform_fresh_search(command)
    else:
        print("An existing record has been found. An incremental blast search will be performed")
        perform_incrermental_search(command, commandOptions, dbinstanceid, current_partsfile)

#incremental_blast("blastn -db ../blastn/whole -query ../blastn/small-query.fasta -outfmt 5 -out ../blastn/result.xml")

def main(search_string):
   incremental_blast(search_string)

if __name__ == "__main__":
    search_string = sys.argv[1]
    print(search_string)
    main(search_string)
