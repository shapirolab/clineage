import sqlite3
from cloud.serialization.cloudpickle import loads
from frogress import bar

SIMULATED_HISTS_PATH = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR3/fastq_human/Calling/ac_mat_2a_prop_100_sim_hists_with_replacement.pickle'

conn = sqlite3.connect('simulations_data.db')
c = conn.cursor()

# Create table
c.execute('''CREATE TABLE simulation
             (id INTEGER PRIMARY KEY AUTOINCREMENT, seed_a integer, seed_b integer, p1 real, cycle integer)''')

c.execute('''CREATE TABLE histogram
             (id INTEGER PRIMARY KEY AUTOINCREMENT, simulation_id INTEGER, xi integer, value real, FOREIGN KEY(simulation_id) REFERENCES simulation(id))''')



f = open(SIMULATED_HISTS_PATH, 'rb').read()
print 'loading existing simulations'
sim_hists = loads(f)

for seeds in bar(sim_hists):
    for seeds_and_proportions in sim_hists[seeds]:
        for cycles in sim_hists[seeds][seeds_and_proportions]:
            gspa = sim_hists[seeds][seeds_and_proportions][cycles]
            seed_a, seed_b = seeds
            c.execute("INSERT INTO simulation VALUES (NULL, {},{},{},{})".format(seed_a, seed_b, seeds_and_proportions[0][1], cycles))
            sim_id = c.lastrowid
            for i in gspa:
                c.execute("INSERT INTO histogram VALUES (NULL, {}, {},{})".format(sim_id, i, gspa[i]))

# Save (commit) the changes
conn.commit()

# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
conn.close()