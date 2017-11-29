import sqlite3
from chemkin import *

db = sqlite3.connect('nucleardb.sqlite')
cursor = db.cursor()

table1 = cursor.execute('''SELECT * FROM ELEMENT_PROPERTIES''').fetchall()
for i in table1:
	print(i)

print('\n')

table2 = cursor.execute('''SELECT * FROM NUCLEAR_EMISSIONS''').fetchall()
for i in table2:
	print(i)


nr = nuclear(['Cf'], ['Pd', 'Te'], [254], [118, 132])
print(nr.check_stable())