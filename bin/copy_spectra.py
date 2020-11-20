#!/usr/bin/env python3

import sys
import sqlite3

target = sys.argv[1]
decoy = sys.argv[2]
setnames = sys.argv[3:]


con = sqlite3.Connection(decoy)

maxmzfrow = con.execute('SELECT MAX(rowid) FROM mzmlfiles').fetchone()[0]
maxmzrow = con.execute('SELECT MAX(rowid) FROM mzml').fetchone()[0]
ioninj, ionmob = '', ''

con.execute('ATTACH ? AS target', (target,))
recs = con.execute('SELECT * FROM target.biosets WHERE set_name IN({})'.format(','.join(['?'] * len(setnames))), tuple(setnames))
con.executemany('INSERT INTO main.biosets VALUES(?, ?)', recs)
recs = con.execute('SELECT * FROM target.mzmlfiles WHERE rowid>?', (maxmzfrow,))
con.executemany('INSERT INTO main.mzmlfiles VALUES(?, ?, ?)', recs)
recs = con.execute('SELECT * FROM target.mzml WHERE rowid>?', (maxmzrow,))
con.executemany('INSERT INTO main.mzml VALUES(?, ?, ?, ?, ?, ?)', recs)
if con.execute('SELECT COUNT(*) FROM target.ioninjtime').fetchone()[0]:
    maxiirow = con.execute('SELECT MAX(rowid) FROM main.ioninjtime').fetchone()[0]
    recs = con.execute('SELECT * FROM target.ioninjtime WHERE rowid>?', (maxiirow,))
    con.executemany('INSERT INTO main.ioninjtime VALUES(?, ?)', recs)
if con.execute('SELECT COUNT(*) FROM target.ionmob').fetchone()[0]:
    maxionrow = con.execute('SELECT MAX(rowid) FROM target.ionmob').fetchone()[0]
    recs = con.execute('SELECT * FROM target.ionmob WHERE rowid>?', (maxionrow,))
    con.executemany('INSERT INTO main.ionmob VALUES(?, ?)', recs)
con.commit()
