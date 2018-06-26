import numpy as np
import os
import pprint

base_path = 'C:\\Users\\matth\\Documents\\realdocs_work\\delft\\schoolwork\\Year 3\\Q4 (DSE)\\final\\FinalVerificationFiles'
fnames = os.listdir(base_path)

output_dic = {}
for fname in fnames:
    dicname = fname[:-4]
    path = os.path.join(base_path, fname)
    some_data = np.loadtxt(path, dtype=str)
    output_dic[dicname] = some_data

pprint.pprint(output_dic)
