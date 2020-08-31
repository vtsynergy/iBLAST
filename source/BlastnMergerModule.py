from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import Bio.Blast.Record 
import math
import numpy as np
import xml.etree.cElementTree as ET
import copy
from timeit import default_timer as timer
import sys

class BlastnMerger:
    def MAX(self, m, n):
        if m >= n:
            return m
        else:
            return n

    def get_record(self, result_file):
        result_handle = open(result_file)
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)
        print("K: ", blast_record.ka_params[1])
        return blast_record


    def compute_length_adjustment(self,
                                 K, 
                                 logK,
                                 alpha_d_lambda,
                                 beta,
                                 query_length,
                                 db_length,
                                 db_num_seqs):

        i = 0                    
        kMaxIterations = 20;     
        m = query_length
        n = db_length
        N = db_num_seqs

        ell = 0.0                
        ss = 0                   
        ell_min = 0
        ell_max = 0 
                    
        converged    = 0
        ell_next = 0 
                     
        length_adjustment = 0
        a  = N
        mb = m * N + n
        c  = n * m - self.MAX(m, n) / K

        if c < 0:
            length_adjustment = 0
            return 1
        else:
            ell_max = 2 * c / (mb + math.sqrt(mb * mb - 4 * a * c))
       
        for i in range(kMaxIterations): 
            ell_bar = 0.0         
            ell      = ell_next
            ss       = (m - ell) * (n - N * ell)
            ell_bar  = alpha_d_lambda * (logK + math.log(ss)) + beta
            if ell_bar >= ell: 
                ell_min = ell
                if ell_bar - ell_min <= 1.0: 
                    converged = 1
                    break
                if ell_min == ell_max: 
                    break  
            else: 
                ell_max = ell
            
            if ell_min <= ell_bar and ell_bar <= ell_max:
                ell_next = ell_bar
            else: 
                ell_next = ell_max if (i == 1)  else (ell_min + ell_max) / 2
            
        if converged:
            length_adjustment = (int) (ell_min)
            ell = math.ceil(ell_min)
            if ell <= ell_max:
              ss = (m - ell) * (n - N * ell)
              if alpha_d_lambda * (logK + math.log(ss)) + beta >= ell:
                length_adjustment = (int) (ell)
        else:
            length_adjustment = (int) (ell_min)
     

        return length_adjustment if converged else 0

    def compute_effective_searchspace(self, db_num_seqs, db_length, query_length):
        alpha = 1.5;
        beta = -2;
        Lambda = 1.28;
        K = 0.46;
        logK = math.log(K);

        alpha_d_lambda = alpha / Lambda;

        length_adjustment = self.compute_length_adjustment(K,
                             logK,
                             alpha_d_lambda,
                             beta,
                             query_length,
                             db_length,
                             db_num_seqs)
        effective_searchspace = (query_length - length_adjustment) * (db_length - db_num_seqs * length_adjustment)
        return effective_searchspace

    def recompute_evalues(self, alignments, effective_searchspace, K, Lambda):
        for alignment in alignments:
            for hsp in alignment.hsps:
                recomputed_evalue = K * effective_searchspace * math.exp(-Lambda * hsp.score)
                if recomputed_evalue < math.pow(10, -180):
                    recomputed_evalue = 0.0
                hsp.expect = recomputed_evalue
        


    def rescale_evalues(self, alignments, part_db_size, full_db_size):
        scaling_factor = (1.0 * full_db_size) / (1.0 * part_db_size)
        print("scaling factor: ", scaling_factor)
        for alignment in alignments:
            for hsp in alignment.hsps:
                rescaled_evalue = hsp.expect * scaling_factor
                hsp.expect = rescaled_evalue
                

    def merge_records(self, record1, record2, K, Lambda, indices):
        scores = []
        indices = []
        #print("Inside merge_records")

        alignments1 = record1.alignments
        alignments2 = record2.alignments

        max_alignment_length = self.MAX(len(alignments1), len(alignments2))
        
        db_size_1 = record1.num_letters_in_database
        db_size_2 = record2.num_letters_in_database
        db_num_seqs_1 = record1.database_sequences
        db_num_seqs_2 = record2.database_sequences

        effective_searchspace = self.compute_effective_searchspace(
            db_num_seqs_1 + db_num_seqs_2, 
            db_size_1 + db_size_2, 
            record1.query_length)
        #print(len(alignments1), len(alignments2))
        if len(alignments1) == 0 and len(alignments2) == 0:
            return alignments1, effective_searchspace, indices;

        if len(alignments1) != 0:
            self.recompute_evalues(alignments1, effective_searchspace, K, Lambda)
        if len(alignments2) != 0:
            self.recompute_evalues(alignments2, effective_searchspace, K, Lambda)

        if len(alignments1) == 0:
            return alignments2, effective_searchspace, indices
        elif len(alignments2) == 0:
            return alignments1, effective_searchspace, indices

        merged_alignments = []

        m = 0
        n = 0

        num_alignments = min(100, len(alignments1) + len(alignments2))
        for k in range(num_alignments):
            if m >= len(alignments1) and n>= len(alignments2):
                break

            min_evalue_1 = 1000000.0
            if m < len(alignments1):
                for hsp in alignments1[m].hsps:
                    if hsp.expect < min_evalue_1:
                        min_evalue_1 = hsp.expect
                        socre_1 = hsp.score
            
            min_evalue_2 = 1000000.0
            if n < len(alignments2):
                for hsp in alignments2[n].hsps:
                    if hsp.expect < min_evalue_2:
                        min_evalue_2 = hsp.expect
                        score_2 = hsp.score

            #print(socre_1, score_2)
            if min_evalue_1 < min_evalue_2:
                merged_alignments.append(alignments1[m])
                scores.append(alignments1[m].hsps[0].score)
                m = m + 1
                indices.append(1)
            elif min_evalue_1 == min_evalue_2:
                if socre_1 >= score_2:
                    merged_alignments.append(alignments1[m])
                    scores.append(alignments1[m].hsps[0].score)
                    m = m + 1
                    indices.append(1)
                else:
                    merged_alignments.append(alignments2[n])
                    scores.append(alignments2[n].hsps[0].score)
                    n = n + 1
                    indices.append(2)
            else:
                merged_alignments.append(alignments2[n])
                scores.append(alignments2[n].hsps[0].score)
                n = n + 1
                indices.append(2)

        #print("Number of alginments chosen from first record: ", m)
        #print("Number of alignments chosen from second record: ", n)
                #print("merged_alignments: ", merged_alignments, "Searchspace: ", effective_searchspace, "Indices: ", indices) 
        #print("alignments: ", merged_alignments)

        #print("searchspace: ", effective_searchspace)

        #print("Indices: ", indices)


        return merged_alignments, effective_searchspace, indices

    def min_evalue_of_alignment(self, alignment):
        min_evalue = 100000.0
        for hsp in alignment.hsps:
            if hsp.expect < min_evalue:
                min_evalue = hsp.expect
        return min_evalue


    def sum_evalue_of_alignment(self, alignment):
        sum_evalue = 0.0
        for hsp in alignment.hsps:
            sum_evalue = sum_evalue + hsp.expect
        return sum_evalue

    def sum_score_of_alignment(self, alignment):
        sum_score = 0.0
        for hsp in alignment.hsps:
            sum_score = sum_score + hsp.score
        return sum_score
                    

    def merge_results(self, current_result_file, delta_result_file):
        current_result_file = "../blastn/current-result.xml"
        delta_result_file = "../blastn/delta-result.xml"
        current_record = self.get_record(current_result_file)
        delta_record = self.get_record(delta_result_file)

        K = current_record.ka_params[1]
        Lambda = current_record.ka_params[0]
        H = current_record.ka_params[2]

        indices = range(500)
        merged_alignments, effective_searchspace = self.merge_records(current_record, delta_record, K, Lambda, indices)
        whole_record = self.get_record("../blastn/whole-result.xml")
        whole_alignments = whole_record.alignments
        effective_searchspace = whole_record.effective_search_space

        print("printing shapes")
        current_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in current_record.alignments])
        delta_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in delta_record.alignments])
        whole_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in whole_record.alignments])
        merged_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in merged_alignments])
        print(merged_evalues.shape)

        print("query length: ", whole_record.query_length)
        for i in range(100):
            current_evalue = self.min_evalue_of_alignment(current_record.alignments[i])
            delta_evalue = self.min_evalue_of_alignment(delta_record.alignments[i])
            merged_evalue = self.min_evalue_of_alignment(merged_alignments[i])
            whole_evalue = self.min_evalue_of_alignment(whole_alignments[i])
            print(current_evalue, ", ", delta_evalue, ", ", merged_evalue, ", ", whole_evalue)

        return current_evalues, delta_evalues, whole_evalues, merged_evalues

    def merge_xml(self, current_xml, delta_xml, merged_alignments, merged_min_evalues):
        merged_hit_list = [ item.hit_id for item in merged_alignments ]
        print("Current file: ", current_xml)
    
        current_tree = ET.parse(current_xml)
        delta_tree = ET.parse(delta_xml)
    
        current_root = current_tree.getroot()
        delta_root = delta_tree.getroot()

        whole_tree = current_tree

        # Remove hits
        current_hit_parent = current_root.find("BlastOutput_iterations").find("Iteration").find("Iteration_hits")
        delta_hit_parent = delta_root.find("BlastOutput_iterations").find("Iteration").find("Iteration_hits")

        current_hits = current_hit_parent.findall("Hit")
        delta_hits = delta_hit_parent.findall("Hit")


        #print(current_hits)
        start = timer()
        for hit in current_hits:
            hit_id = hit[1].text
            #print(hit_id)
            if hit_id not in merged_hit_list:
                current_hit_parent.remove(hit)
                
                
        for hit in delta_hits:
            hit_id = hit[1].text
            if hit_id not in merged_hit_list:
                delta_hit_parent.remove(hit)
        end = timer()
        print("Elapsed Time: ", end - start)

        current_hits = current_hit_parent.findall("Hit")
        delta_hits = delta_hit_parent.findall("Hit")

        print(len(delta_hits))
        current_hit_ids = [ hit[1].text for hit in current_hits ]
        delta_hit_ids = [ hit[1].text for hit in delta_hits ]

        print("current: ", len(current_hit_ids))
        print("Delta: ", len(delta_hit_ids))
        print("Merged: ", len(merged_hit_list))

        c = 0
        d = 0
        m = 0
        start = timer()
        for m in range(len(merged_hit_list)):
            print(m, len(current_hits))
            hit_id = merged_hit_list[m]
            print("d: " + str(d))
            if len(delta_hit_ids) != 0:
                if hit_id == delta_hit_ids[d]:
                    delta_hits[d].find("Hit_num").text = str(m + 1)
                    current_hit_parent.insert(m, delta_hits[d])
                    current_hits.insert(m, delta_hits[d])
                    current_hit_ids.insert(m, delta_hits[d][1].text)
                    d += 1

            if hit_id == current_hit_ids[m]:
                current_hits[m].find("Hit_num").text = str(m + 1)

        for i in range(len(current_hits)):
            #print(current_hits[i].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text)
            current_hits[i].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text = str(merged_min_evalues[i])

        whole_tree.write("new_merged_result.xml")

        end = timer()
        print("Elapsed Time: ", end - start)    

        print("Elapsed Time: ", end - start)

    def merge(self, current_result_file, delta_result_file):
        current_record = self.get_record(current_result_file)
        delta_record = self.get_record(delta_result_file)

        K = current_record.ka_params[1]
        Lambda = current_record.ka_params[0]
        H = current_record.ka_params[2]

        indices = list(range(500))
        merged_alignments, effective_searchspace, indices = self.merge_records(current_record, delta_record, K, Lambda, indices)

        current_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in current_record.alignments])
        delta_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in delta_record.alignments])
        merged_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in merged_alignments])
        print(merged_evalues.shape)

        for i in range(100):
            current_evalue = self.min_evalue_of_alignment(current_record.alignments[i])
            delta_evalue = self.min_evalue_of_alignment(delta_record.alignments[i])
            merged_evalue = self.min_evalue_of_alignment(merged_alignments[i])
            print(current_evalue, ", ", delta_evalue, ", ", merged_evalue)


        self.merge_xml(current_result_file, delta_result_file, merged_alignments, merged_evalues)
        return current_evalues, delta_evalues, merged_evalues

    def get_all_records(self, result_file):
        blast_records = []
        result_handle = open(result_file)
        blast_record_iter = NCBIXML.parse(result_handle)


        while(True):
            try:
                blast_record = next(blast_record_iter)
                blast_records.append(blast_record)
            except:
                break
        result_handle.close()
        print("Number of records: " + repr(len(blast_records)))
        return blast_records



    def merge_multi(self, current_result_file, delta_result_file, output_file):
        print("Received output file name: ", output_file)
        current_records = self.get_all_records(current_result_file)
        #print(current_records)
        delta_records = self.get_all_records(delta_result_file)

        merged_records = []
        #print("PRINTING KA_PARAMS", current_records[0].ka_params)
        K = current_records[0].ka_params[1]
        Lambda = current_records[0].ka_params[0]
        H = current_records[0].ka_params[2]

        indices = list(range(500))

        current_query_to_record = {}
        delta_query_to_record = {}

        for record in current_records:
            current_query_to_record[record.query] = record
        for record in delta_records:
            delta_query_to_record[record.query] = record

        num_queries = len(current_records)

        current_tree = ET.parse(current_result_file)
        delta_tree = ET.parse(delta_result_file)
    
        current_root = current_tree.getroot()
        delta_root = delta_tree.getroot()

        current_iterations = current_root.find("BlastOutput_iterations").findall("Iteration")
        delta_iterations = delta_root.find("BlastOutput_iterations").findall("Iteration")

        #print("Iteration lengths: ", len(current_iterations), len(delta_iterations))

        
        # Set up output XML file
        merged_root = ET.Element("BlastOutput")
        for i in range(8):
            merged_root.insert(i, copy.deepcopy(current_root[i]))
        merged_iterations_root = ET.SubElement(merged_root, "BlastOutput_iterations")

        merged_iterations = []

        current_statistics = current_iterations[0].find("Iteration_stat").find("Statistics")
        delta_statistics = delta_iterations[0].find("Iteration_stat").find("Statistics")

        db_num_c = int(current_statistics.find("Statistics_db-num").text) 
        db_num_d = int(delta_statistics.find("Statistics_db-num").text)
        db_num_m = db_num_c + db_num_d

        db_len_c = int(current_statistics.find("Statistics_db-len").text) 
        db_len_d = int(delta_statistics.find("Statistics_db-len").text)
        db_len_m = db_len_c + db_len_d

        for q in range(num_queries):
            #print("Processing query# ", q)
            current_record = current_records[q]
            delta_record = delta_records[q]
            merged_alignments, effective_searchspace_m, merged_indices = self.merge_records(current_record, delta_record, K, Lambda, indices)
            #print("Number of merged alignments: " + repr(len(merged_alignments)))
            
            #print(q, len(current_record.alignments), len(delta_record.alignments), len(merged_alignments))
            #current_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in current_record.alignments])
            #delta_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in delta_record.alignments])
            merged_evalues = np.asarray([self.min_evalue_of_alignment(alignment) for alignment in merged_alignments])
            #print(merged_evalues)
            current_iteration = current_iterations[q]
            delta_iteration = delta_iterations[q]    

            merged_iteration = self.merge_two_iterations(current_iteration, delta_iteration, merged_alignments, merged_evalues, merged_indices)
            merged_iterations_root.insert(q, merged_iteration)
            merged_iterations.append(merged_iteration)

            merged_statistics = merged_iteration.find("Iteration_stat").find("Statistics")

            merged_statistics.find("Statistics_db-num").text = str(db_num_m)
            merged_statistics.find("Statistics_db-len").text = str(db_len_m)
            merged_statistics.find("Statistics_eff-space").text = str(effective_searchspace_m)

        print("Writing to Output File: ", output_file)
        ET.ElementTree(merged_root).write(output_file, encoding="utf-8", xml_declaration=True)
        print("Done writing. going back to entry point.")
        #ET.ElementTree(merged_root).write("../blastn/failsafe.xml")
        

    def merge_two_iterations(self, current_iteration, delta_iteration, merged_alignments, merged_min_evalues, merged_indices):
        
        scores = copy.deepcopy([alignment.hsps[0].score for alignment in merged_alignments])

        merged_hit_ids = [ item.hit_id for item in merged_alignments ]
        #print("Length of merged hit ids: ", len(merged_hit_ids))
    
        merged_iteration = ET.Element("Iteration")
        ci = 0
        #print("Printing children")
        for child in current_iteration:
            if child.tag != "Iteration_hits":
                #print(ci, child.tag)
                merged_iteration.insert(ci, child)
            ci += 1

        merged_iteration_hits = ET.Element("Iteration_hits")
        #merged_iteration.insert(3, merged_iteration_hits)

    
        current_hits = list(current_iteration.find("Iteration_hits").findall("Hit"))
        delta_hits = list(delta_iteration.find("Iteration_hits").findall("Hit"))
        #print("Current, delta, and merged hit lengths:", len(current_hits), len(delta_hits), len(merged_alignments))
        #print("Current hits: ", current_hits)
        #print("Delta hits: ", delta_hits)

        c = 0
        d = 0

        #print("Printing scores again")
        #print(scores)

        for m in range(len(merged_alignments)):
            #print(repr(m), repr(scores[m]))


            if c < len(current_hits) and current_hits[c].find("Hit_id").text == merged_hit_ids[m]:
                merged_hit = copy.deepcopy(current_hits[c])
                merged_iteration_hits.insert(m, merged_hit)
                #score = current_hits[c].find("Hit_hsps").find("Hsp").find("Hsp_score").text
                merged_iteration_hits[m].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text = str(merged_min_evalues[m]) 
                #print(repr(m) + " : " + repr(merged_alignments[m].hsps[0].score) , " <> ", score, " <> ", scores[m])
                c += 1

                #print(1, merged_indices[m], merged_min_evalues[m], merged_iteration_hits[m].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text)

            elif d < len(delta_hits) and delta_hits[d].find("Hit_id").text == merged_hit_ids[m]:
                merged_hit = copy.deepcopy(delta_hits[d])
                merged_iteration_hits.insert(m, merged_hit)
                #score = delta_hits[d].find("Hit_hsps").find("Hsp").find("Hsp_score").text
                merged_iteration_hits[m].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text = str(merged_min_evalues[m]) 
                #print(repr(merged_alignments[m].hsps[0].score) , " <> ", score, " <> ", scores[m])
                #print(2, merged_indices[m], merged_min_evalues[m], merged_iteration_hits[m].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text)
                d += 1
            else:
                print("Nothing happened")

        #for i in range(500):
        #    print(merged_iteration_hits[i].find("Hit_hsps").find("Hsp").find("Hsp_evalue").text, merged_min_evalues[i])

        merged_iteration.insert(4, merged_iteration_hits)


        
        #print("Printing merged min evalues:")
        #for meval in merged_min_evalues:
        #        print(meval)
        return merged_iteration


def main(argv):
    result1 = argv[1]
    result2 = argv[2]
    combined_result = argv[3]
    merger = BlastnMerger()
    merger.merge_multi(result1, result2, combined_result)
    

if __name__ == "__main__":
    main(sys.argv)

#merger.merge_multi("../blastn/current-30-result-500q.xml", "../blastn/delta-30-40-result-500q.xml", "../blastn/derived-whole-40-result-500q-new.xml")

#python BlastnMergerModule.py "../blastn/current-30-result-500q.xml" "../blastn/delta-30-40-result-500q.xml" "../blastn/derived-whole-40-result-500q-new.xml"










