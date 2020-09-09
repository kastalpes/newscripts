from GL import *
from collections import defaultdict

def return_average_quantities(file_list=[], all_quan = defaultdict(lambda: list())):
    for this_name in file_list:
        if not os.path.exists(this_name):
            print("No such file, skipping: %s"%this_name)
            continue
        fptr = h5py.File(this_name,'r')
        for key in fptr:
            mything =  list(fptr[key][:].flatten())
            all_quan[key] += mything
        try:
            pass
        except:
            raise
        finally:
            fptr.close()
    for key in all_quan:
        all_quan[key] = np.array( all_quan[key]).flatten()
    return all_quan

def files_from_output(path_to_output_log, suffix=".AverageQuantities.h5"):
    if 'OutputLog' not in path_to_output_log.split("/"):
        output_log = path_to_output_log + "/OutputLog"
        path = path_to_output_log
    else:
        output_log=path_to_output_log
        path = "/".join(output_log.split("/")[:-1])

    fptr = open(output_log,'r')
    lines = fptr.readlines()
    fptr.close()
    files = []
    for line in lines:
        files.append( path+line.split()[2][1:]+suffix)
    return files

def all_quan_from_outputlog(output_log):
    file_list= files_from_output(output_log)
    

    quan = return_average_quantities(file_list)
    return quan

def all_quan_from_taxi(car):
    dumb=[]
    file_list=[]
    for n in car.frame_dict:
        file_list.append("%s/%s.AverageQuantities.h5"%(car.directory,car.frame_dict[n]['dsname']))
    all_quan = return_average_quantities(file_list)
    return all_quan

def all_quan_from_files(direc,frame_list):
    file_list=[]
    for frame in frame_list:
        file_list.append("%s/data%04d.AverageQuantities.h5"%(direc,frame))
        print("%s/data%04d.AverageQuantities.h5"%(direc,frame))
    print(file_list)
    pdb.set_trace()
    all_quan = return_average_quantities(file_list)
    return all_quan
