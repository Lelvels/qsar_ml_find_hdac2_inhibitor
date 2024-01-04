import pandas as pd
import numpy as np
from tabulate import tabulate

import os
import psycopg2

from dotenv import load_dotenv
import logging
from tqdm import tqdm
import multiprocessing
import time

#Rdkit ultis
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import MACCSkeys


"""Notification: This file is used to calculate the average distance between a smiles and its 5 nearest-neighbor in the database."""

#Logging init
logging.basicConfig(filename='../../output/logs/hdac2-average-distance-auto-update.txt', level=logging.INFO)

#Ultis
class PubchemScreeningDAO:    
    """
    Data access object of the data after preprocessed the PubChem database, 
        the table include 4 columns: id, cid, smiles, molecular_weight
    """  
    CREATE_PREPROCESSING_TABLE = """CREATE TABLE IF NOT EXISTS pubchem_compound_preprocessed (
        id SERIAL PRIMARY KEY,
        cid INTEGER UNIQUE NOT NULL, 
        smiles TEXT UNIQUE NOT NULL,
        molecular_weight decimal NOT NULL,
        average_distance float(24) NOT NULL
    );
    CREATE INDEX IF NOT EXISTS idx_cid ON pubchem_compound_preprocessed (cid);
    """
    SELECT_DATA_BETWEEN_CID = """SELECT * FROM 
        pubchem_compound_preprocessed WHERE molecular_weight > 200 and cid >= %s and cid <= %s;"""
    
    INSERT_COMPOUND = "INSERT INTO pubchem_compound_preprocessed (cid, smiles, molecular_weight, average_distance) VALUES (%s, %s, %s, %s);"

    SELECT_FIRST_DATA = "SELECT * FROM pubchem_compound_preprocessed LIMIT %s;"

    SELECT_DATA_BY_CID = """SELECT * FROM pubchem_compound_preprocessed WHERE cid = %s;"""

    SELECT_NUMBERS_OF_ROWS = """SELECT COUNT(*) as row_count FROM pubchem_compound_preprocessed;"""

    MAX_CID_QUERY = """SELECT MAX(cid) FROM pubchem_compound_preprocessed;"""

    MIN_CID_QUERY = """SELECT MIN(cid) FROM pubchem_compound_preprocessed;"""

    DELETE_COMPOUND = """DELETE FROM pubchem_compound_preprocessed WHERE cid = %s;"""

    UPDATE_COMPOUND = """UPDATE pubchem_compound_preprocessed \
        SET smiles=%s, molecular_weight=%s, average_distance=%s
        WHERE cid = %s;"""
    
    UPDATE_AVG_DIST_COMPOUND = """UPDATE pubchem_compound_preprocessed \
        SET average_distance=%s
        WHERE cid = %s;"""
    
    MAX_CID_QUERY = "SELECT MAX(cid) FROM pubchem_compound_preprocessed"

    MIN_CID_QUERY = "SELECT MIN(cid) FROM pubchem_compound_preprocessed"
    connection = None
    max_cid = None
    min_cid = None
    
    def __init__(self):
        load_dotenv("../env/.env")
        self.connection = psycopg2.connect(os.environ['DATABASE_URL'])
        self.get_min_max_cid()
        logging.info(f"Min CID - Max CID: {self.min_cid} - {self.max_cid}")

    def get_min_max_cid(self):
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(self.MIN_CID_QUERY)
                self.min_cid = cursor.fetchall()[0][0]
                cursor.execute(self.MAX_CID_QUERY)
                self.max_cid = cursor.fetchall()[0][0]
                return self.min_cid, self.max_cid

    def create_tables(self):
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(self.CREATE_PREPROCESSING_TABLE)
    
    def read_first_rows(self, number_of_data):
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(self.SELECT_FIRST_1000_ALL_DATA, number_of_data)
                column_names = [desc[0] for desc in cursor.description]
                data_list = cursor.fetchall()
                return pd.DataFrame(data_list, columns=column_names)
    
    def get_data_between_cid(self, first_cid, second_cid):
        if(first_cid < self.min_cid or second_cid > self.max_cid):
            logging.info("Invalid cid!")
            return
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(self.SELECT_DATA_BETWEEN_CID, (first_cid, second_cid))
                data_list = cursor.fetchall()
                column_names = [desc[0] for desc in cursor.description]
                return pd.DataFrame(data_list, columns=column_names)
    
    def insert(self, cid, smiles, molecular_weight, average_distance):
        with self.connection:
            with self.connection.cursor() as cursor:
                try:
                    cursor.execute(self.INSERT_COMPOUND, (cid, smiles, molecular_weight, average_distance))
                    return True
                except psycopg2.IntegrityError as error:
                    if 'duplicate key value violates unique constraint' in str(error):
                        logging.error(f"Skipping data with cid={cid} due to unique constraint violation.")
                    else:
                        logging.error(str(error))
                        raise  # Re-raise the exception if it's not related to unique constraint violation
                except Exception as e:
                    logging.error("An exception occurred: " + str(e))
        return False  # Return False if no successful insertion or if an exception occurs    

    def update(self, cid, smiles, molecular_weight, average_distance):
        with self.connection:
            with self.connection.cursor() as cursor:
                try:
                    cursor.execute(self.UPDATE_COMPOUND, (int(cid), str(smiles), float(molecular_weight), float(average_distance)))
                    return True
                except psycopg2.Error as error:
                    if 'duplicate key value violates unique constraint' in str(error):
                        logging.error(f"Skipping data with cid={cid} due to unique constraint violation.")
                    else:
                        logging.error(str(error))
                        raise
                except Exception as e:
                    logging.error("An exception occurred: " + str(e))
        return False  # Return False if no successful insertion or if an exception occurs
    
    def update(self, cid, smiles, molecular_weight, average_distance):
        try:
            cid = int(cid)
            smiles = str(smiles)
            molecular_weight = float(molecular_weight)
            average_distance = float(average_distance)
            with self.connection:
                with self.connection.cursor() as cursor:
                    cursor.execute(
                        self.UPDATE_COMPOUND,
                        (smiles, molecular_weight, average_distance, cid)  # Add the CID as the last parameter for WHERE clause
                    )
            return True
        except psycopg2.Error as error:
            if 'duplicate key value violates unique constraint' in str(error):
                logging.error(f"Skipping data with cid={cid} due to unique constraint violation.")
            else:
                logging.error(str(error))
            return False
        except Exception as e:
            logging.error("An exception occurred: " + str(e))
        return False
    
    def update_average_distance(self, cid, average_distance):
        average_distance = float(average_distance)
        cid = int(cid)
        with self.connection:
            with self.connection.cursor() as cursor:
                try:
                    cursor.execute(self.UPDATE_AVG_DIST_COMPOUND, (average_distance, cid))
                    return True
                except psycopg2.Error as error:
                    if 'duplicate key value violates unique constraint' in str(error):
                        logging.error(f"Skipping data with cid={cid} due to unique constraint violation.")
                    else:
                        logging.error(str(error))
                        raise
                except Exception as e:
                    logging.error("An exception occurred: " + str(e))
        return False
    
    def delete(self, cid):
        with self.connection:
            with self.connection.cursor() as cursor:
                try:
                    cursor.execute(self.DELETE_COMPOUND, (cid))
                    return True
                except psycopg2.Error as error:
                    logging.error(str(error))
                except Exception as e:
                    logging.error("An exception occurred: " + str(e))
        return False
    
    def get_number_of_rows(self):
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(self.SELECT_NUMBERS_OF_ROWS)
                nor = cursor.fetchall()
                return nor
    
    def create_connection(self):
        load_dotenv("../env/.env")
        self.connection = psycopg2.connect(os.environ['DATABASE_URL'])

    def close_connection(self):
        self.connection.close()

#Morgan2
def morgan_fpts(data):
    Morgan_fpts = []
    count = 0
    with tqdm(total=len(data), desc='Morgan2 encoding:') as pbar:
        for i in data:
            try:
                mol = Chem.MolFromSmiles(i)
            except:
                print("An exception occurred with " + str(count))
                continue
            fpts = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
            mfpts = np.array(fpts)
            Morgan_fpts.append(mfpts)
            count += 1
            pbar.update(1)  # Update the progress bar
    return np.array(Morgan_fpts)

def import_train_fpts():
    train_test_path = "../../data_for_modeling/train_test_data/v3/HDAC2_unclean_data_logistic_regression.xlsx"
    train_dataset = pd.read_excel(train_test_path, sheet_name='train_dataset')
    print("Train data imported:")
    print(len(train_dataset))
    X_Train = morgan_fpts(train_dataset['SMILES'])
    return X_Train

def cal_tc(array1, array2):
    if len(array1) != len(array2):
        raise ValueError("The arrays must have the same length.")
    
    # Calculate the Tanimoto coefficient
    intersection = sum(a and b for a, b in zip(array1, array2))
    union = sum(a or b for a, b in zip(array1, array2))
    
    if union == 0:  # Handle the case when both arrays are all zeros
        return 0.0
    else:
        tanimoto_coefficient = intersection / union
        return tanimoto_coefficient

def find_nearest_neighbors_distance_single_core(X_Screening, n_neighbors, X_Train):
    nearest_neighbors_distances = []
    nearest_neighbors_indices = []
    X_Screening = np.array(X_Screening)
    
    if X_Screening.shape[1] != X_Train.shape[1]:
        raise ValueError("X_Screening bit vectors must have the same size as X_Train bit vectors: " + str(X_Train.shape[1]))
    
    for screening_vector in X_Screening:
        tc_dist = []
        for X_curr in X_Train:
            tc_dist.append(cal_tc(screening_vector, X_curr))
        
        # Get the indices of the first n_neighbors elements with the largest Tanimoto coefficients
        nearest_neighbor_indices = sorted(range(len(tc_dist)), key=lambda i: tc_dist[i], reverse=True)[:n_neighbors]
        
        # Extract the distances to the nearest neighbors
        nearest_neighbors_distances.append([tc_dist[i] for i in nearest_neighbor_indices])
        nearest_neighbors_indices.append(nearest_neighbor_indices)
    
    return nearest_neighbors_distances, nearest_neighbors_indices  

def calculate_distances(identifier, screening_vectors, X_Train):
    distances = []
    for screening_vector in screening_vectors:
        tc_dist = np.sum(screening_vector & X_Train, axis=1) / np.sum(screening_vector | X_Train, axis=1)
        distances.append(tc_dist)
    return identifier, distances

def find_nearest_neighbors_distance(X_Screening, n_neighbors, X_Train):
    nearest_neighbors_distances = []
    nearest_neighbors_indices = []

    X_Screening = np.array(X_Screening)
    X_Train = np.array(X_Train)
    if(np.size(X_Screening) == 0):
        return nearest_neighbors_distances, nearest_neighbors_indices
     
    if X_Screening.shape[1] != X_Train.shape[1]:
        raise ValueError("X_Screening bit vectors must have the same size as X_Train bit vectors: " + str(X_Train.shape[1]))

    num_processes = multiprocessing.cpu_count()
    
    if len(X_Screening) <= num_processes:
        screening_chunks = [(i, X_Screening[i:i + 1]) for i in range(len(X_Screening))]
    else:
        chunk_size = len(X_Screening) // num_processes
        screening_chunks = [(i, X_Screening[i:i + chunk_size]) for i in range(0, len(X_Screening), chunk_size)]

    pool = multiprocessing.Pool(processes=num_processes)
    results = pool.starmap(calculate_distances, [(i, chunk, X_Train) for i, chunk in screening_chunks])
    pool.close()
    pool.join()

    # Sort the results by identifier to ensure the correct order
    results.sort(key=lambda x: x[0])

    # Extract the distances and indices
    for _, distances in results:
        for distance in distances:
            # Get the indices of the first n_neighbors elements with the largest Tanimoto coefficients
            nearest_neighbor_indices = np.argsort(distance)[::-1][:n_neighbors]

            # Extract the distances to the nearest neighbors
            nearest_neighbors_distances.append([distance[i] for i in nearest_neighbor_indices])
            nearest_neighbors_indices.append(nearest_neighbor_indices)

    return nearest_neighbors_distances, nearest_neighbors_indices

def update_average_distance(screening_dao, n_neighbors, start_cid, end_cid, X_Train):
    logging.info(f"[+] Update average distance for: {start_cid} to {end_cid}")
    working_dataset = screening_dao.get_data_between_cid(start_cid, end_cid)
    if(len(working_dataset) > 0):
        X_Screening = morgan_fpts(working_dataset['smiles'])
        logging.info(f"[-] Start finding nearest neighbor for: {len(X_Screening)}")
        #Find nearest neighbor
        dist_array, nn_idx = find_nearest_neighbors_distance(X_Screening=X_Screening, n_neighbors=n_neighbors, X_Train=X_Train)
        #Inserted counts
        insert_counts = 0
        logging.info("[-] Finished finding, update in database!")
        for idx, row in working_dataset.iterrows():
            cid = row['cid']
            nn_dist = dist_array[idx]
            #Calculate average distance
            avg_distance = np.average(nn_dist)
            #Logging data
            # logging.info("Nearest neighbor index with cid: " + str(row['cid']) + ": " + str(nn_idx[idx]))
            # logging.info("Average distance: " + str(avg_distance))        
            #Update in database
            result = screening_dao.update_average_distance(cid=cid, average_distance=avg_distance)
            insert_counts = insert_counts + result
        logging.info("[-] Update in database: " + str(insert_counts))
        logging.info("\n")
    else:
        logging.info("[-] Empty data, skip this batch!")
        print("[-] Empty data, skip this batch!")

def main():
    #Import the X_Train data
    X_Train = import_train_fpts()
    #Init database dao
    screening_dao = PubchemScreeningDAO()
    #Update the average data
    # Choosing the starting CID to process by min_cid, and the final cid by max_cid.
    # Steps: the number of precessed data for each batch
    # nn_nums: numbers of nearest-neighbors
    nn_nums = 5
    min_cid = 140 * 10 ** 6
    max_cid = 140 * 10 ** 6+1
    step = 1
    
    cid_range = range(int(min_cid), int(max_cid), step)
    for i in cid_range:
        start_cid = i
        end_cid = i + step - 1
        # Starting
        logging.info("\n")
        logging.info(f"[+] Starting new update, from {start_cid} to {end_cid}!")
        print(f"[+] Starting new update, from {start_cid} to {end_cid}!")
        update_average_distance(screening_dao=screening_dao, n_neighbors=nn_nums, start_cid=start_cid, end_cid=end_cid, X_Train=X_Train)
        print("[+] Finished, sleep 5s to regulate then continue!")
        time.sleep(5)

    #Close database connection
    screening_dao.close_connection()

if __name__=="__main__":
    main()