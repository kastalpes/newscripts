#!/usr/bin/env python
"""Scrubs the input file for parameter filename, and [initial,final]\times[time,cycle,dt]."""

import os
#exec(compile(open('%s/yt3_scripts/go_lite'%os.environ['HOME']).read(), '%s/yt3_scripts/go_lite'%os.environ['HOME'], 'exec'))
from optparse import OptionParser
import glob
import sys
from GL import *
import re
import pdb
import numpy as np
nar = np.array
from datetime import *
from xml.dom import minidom


try:
    import cadac
    from . import RunTrackerStuff
    sysInfo = RunTrackerStuff.sysInfo() #Get info about the system and user.
    sys.path.append(sysInfo.cadacInstall)
    RunTracker_on = True
except:
    RunTracker_on = False
    sysInfo = None

class stepInfo:
    dt = -1
    time = -1
    cycle = -1
    wall = None
    set = False
parser = OptionParser("ParseTimesteps output")
parser.add_option("-o","--stdout", dest="stdout", action="store",default=None,type="string")
parser.add_option("-e","--stderr", dest="stderr", action="store",default=None,type="string")
parser.add_option("-j","--jobID", dest="jobID", action="store",default=None,type="string")
(options, args) = parser.parse_args()

#I started trying to discriminate between stderr and stdout,
#but I don't think it matters.
if options.stderr:
    if glob.glob(options.stderr) ==[]:
        sys.stderr.write("ParseTimesteps.py file "+options.stderr+" not found.\n")
        sys.exit(1)
    filename = options.stderr


def scrubParameterFile(string):
    """Matches each line in the file *string* (good name, eh?) to the regular expressions
    in re_list.
    This isn't presently used."""
    
    if string == None:
        return None
    print("Param File To Open: ", string)
    if glob.glob(string) == []:
        return None

    fptr = open(string, 'r')
    pa = []
    re_list = []
    re_list.append(re.compile(r'\s*(TopGridDimensions)\s*=\s*(\d*)\s*(\d*)\s*(\d*)'))

    for line in fptr:
        if len(re_list) == 0:
            break
        for reg in re_list:
            match = reg.match(line)
            if match != None:
                re_list.remove(reg) #save us some effort.
                #This is where you should put the new parameters.
                if match.group(1) == 'TopGridDimensions':
                    pa.append(cadac.Parameter('TopGridDimensions0',match.group(2)))
                    if match.group(3) != '':
                        pa.append(cadac.Parameter('TopGridDimensions1',match.group(3)))
                    if match.group(4) != '':
                        pa.append(cadac.Parameter('TopGridDimensions3',match.group(4)))
        
                elif match.group(1) == "monkey":
                    monkey = "wtf?"
    plist = cadac.ParameterList(pa)
    fptr.close()
    return plist

def scrub_log_file(filename, all_output=None):
    counter_xx  = 0
    counter_xxx = 0
    fptr = open(filename,'r')

    initialFile_start = re.compile(r'.*Successfully read in parameter file\s+(\S+)\.+')
    initialFile_restart = re.compile(r'.*Successfully read ParameterFile\s+(\S+).+')
    initialFileString = None

    oldTicker = re.compile(r'.*TopGrid\s*cycle\s*=\s*(\S+)\s*dt\s*=\s*(\S*)\s*time\s*=\s*(\S+)\s*(wall = (\S+))?')
    newTicker = re.compile(r'.*STEP_INFO\s*N\s*=\s*(\S+)\s*\s*DT\s*=\s*(\S+)\s*,?\s*T\s*=\s*(\S+)')
    weekOfCodeTicker_with_wall = re.compile(r'.*TopGrid\s*dt\s*=\s(\S+)\s*time\s*=\s*(\S+)\s*cycle\s*=\s*(\S+)\s*(wall = (\S+))?')
    #weekOfCodeTicker = re.compile(r'.*TopGrid\s*dt\s*=\s(\S+)\s*time\s*=\s*(\S+)\s*cycle\s*=\s*(\S+)')
    weekOfCodeTicker = re.compile(r'.*TopGrid\s*dt\s*=\s(\S+)\s*time\s*=\s*(\S+)\s*cycle\s*=\s*(\S+)\s*z\s*=\s*(\S+)')
#weekOfCodeTicker = re.compile(r'.*TopGrid\s*dt\s*=\s(\S+)\s*time\s*=\s*(\S+)\s*cycle\s*=\s*(\S+)\s*(wall = (\d*\.\d*))?C.*)?')
#TopGrid dt = 1.690168e-06     time = 0.025478588    cycle = 10500
    old_map = {'cycle':1,'dt':2,'time':3,'wall':5}
    woc_map_wall = {'cycle':3,'dt':1,'time':2,'wall':5}
    woc_map = {'cycle':3,'dt':1,'time':2,'wall':5}
    TickerList = [(woc_map,weekOfCodeTicker), (old_map,oldTicker), (old_map,newTicker),(woc_map_wall, weekOfCodeTicker_with_wall)]
    if all_output is None:
        initialStepInfo = stepInfo()
        finalStepInfo = stepInfo()
        firstTime=True
        all_steps = {'cycle':[], 'dt':[], 'time':[],'output':[],'input':[]}
    else:
        initialStepInfo = all_stuff['initialStepInfo']
        finalStepInfo = all_stuff['finalStepInfo']
        firstTime=False
        all_steps = all_stuff['all_steps']
    last_ten_lines = []
    for line in fptr: 

        last_ten_lines.append(line)
        if len(last_ten_lines) == 11:
            last_ten_lines.pop(0)
        match = initialFile_start.match(line)
        initialize_line=False
        if match != None:
            initialFileString = match.group(1)
            initialize_line=True
        
        match = initialFile_restart.match(line)
        if match != None:
            initialFileString = match.group(1)
            initialize_line=True
        if initialize_line:
            print(fname, initialFileString)
            re_input = re.compile(r'([^\d]+)(\d\d\d\d)/([^\d]*)')
            match_restart = re_input.match(initialFileString)
            if match_restart is not None:
                read_directory=match_restart.group(1)
                read_number   =int(match_restart.group(2))
                read_fname    =match_restart.group(3)
                this_time = -1
                this_cycle = -1
                if len(all_steps['cycle']) > 0:
                    this_time = all_steps['time'][-1]
                    this_cycle = all_steps['cycle'][-1]
                input_dict={'dir':read_directory,'number':int(read_number),'fname':read_fname,
                             'cycle':int(this_cycle),'time':float(this_time)}
                all_steps['input'].append(input_dict)
                #print "NNNNNNNNNNN",input_dict 
            #print "WOAH", line, initialFileString, [val['number'] for val in all_steps['input']]
            #print input_dict
            #print "WOOOOH"
            continue


        
        for reg_map,ticker in TickerList:
            match = ticker.match(line)
            if match != None:
                if firstTime == True:
                    initialStepInfo.cycle = int(match.group(reg_map['cycle']))
                    initialStepInfo.dt =    float(match.group(reg_map['dt']))
                    initialStepInfo.time =  float(match.group(reg_map['time']))
                    initialStepInfo.set = True
                    firstTime = False
                    all_steps['cycle'].append(initialStepInfo.cycle)
                    all_steps['dt'].append(initialStepInfo.dt)
                    all_steps['time'].append(initialStepInfo.time)
                    try:
                        initialStepInfo.wall = float(match.group(reg_map['wall']))
                    except:
                        pass
                    TickerList = [(reg_map,ticker)]
                else:
                    finalStepInfo.cycle =int( match.group(reg_map['cycle']))
                    finalStepInfo.dt =   float( match.group(reg_map['dt']))
                    finalStepInfo.time = float( match.group(reg_map['time']))
                    finalStepInfo.set = True
                    all_steps['cycle'].append(finalStepInfo.cycle)
                    all_steps['dt'].append(finalStepInfo.dt)
                    all_steps['time'].append(finalStepInfo.time)
                    try:
                        finalStepInfo.wall = float(match.group(reg_map['wall']))
                    except:
                        pass
        re_data = re.compile(r'DATA dump: ./([^\d]+)(\d\d\d\d)/([^\d]*)')
        match_DD = re_data.match(line)
        if match_DD is not None:
            dump_directory=match_DD.group(1)
            dump_number   =int(match_DD.group(2))
            dump_fname    =match_DD.group(3)
            this_time = -1
            this_cycle = -1
            if len(all_steps['cycle']) > 0:
                this_time = all_steps['time'][-1]
                this_cycle = all_steps['cycle'][-1]
            output_dict={'dir':dump_directory,'number':int(dump_number),'fname':dump_fname,
                         'cycle':int(this_cycle),'time':float(this_time)}
            all_steps['output'].append(output_dict)

    counter_xxx += 1
    fptr.close()



    if RunTracker_on:
        ua=[]
        if initialFileString != None:
            ua.append(cadac.UserField('startFile',initialFileString))
        if initialStepInfo.set and not finalStepInfo.set:
            #then the first step wasn't finished.
            finalStepInfo = initialStepInfo
        ua.append(cadac.UserField('init_time',initialStepInfo.time ))
        ua.append(cadac.UserField('final_time',finalStepInfo.time ))

        ua.append(cadac.UserField('init_cycle',initialStepInfo.cycle))
        ua.append(cadac.UserField('final_cycle',finalStepInfo.cycle ))

        ua.append(cadac.UserField('init_dt',initialStepInfo.dt ))
        ua.append(cadac.UserField('final_dt',finalStepInfo.dt ))

        userlist = cadac.UserFieldList(ua)

        parameterList = None #scrubParameterFile(initialFileString)

        output_dom = minidom.parseString('<appendMe></appendMe>')

        list_dom= minidom.parseString(userlist.toxml()).firstChild
        output_dom.firstChild.appendChild(list_dom)

        if parameterList != None:
            param_dom = minidom.parseString(parameterList.toxml()).firstChild
            output_dom.firstChild.appendChild(param_dom)

        #Sanatize wall info. This should be in the regexp, but my regexp-fu is weak today.
        if 'C' in finalStepInfo.wall:
            finalStepInfo.wall = finalStepInfo.wall[ :finalStepInfo.wall.index('C') ]
        if 'C' in initialStepInfo.wall:
            initialStepInfo.wall = initialStepInfo.wall[ :initialStepInfo.wall.index('C') ]

        if initialStepInfo.wall != None:
            #This sure is a lot of work for one line of xml.
            #there has to be an easier way...
            round = int(float(initialStepInfo.wall))  #really don't care about microseconds.

            iso_string = datetime.fromtimestamp(round).isoformat()
            start_xml=minidom.parseString('<StartTime>'+iso_string+'</StartTime>').firstChild
            output_dom.firstChild.appendChild( output_dom.createTextNode('\n'))
            output_dom.firstChild.appendChild(start_xml)

        if finalStepInfo.wall != None:
            round = int(float(finalStepInfo.wall))  #really don't care about microseconds.
            iso_string = datetime.fromtimestamp(round).isoformat()
            end_xml=minidom.parseString('<EndTime>'+iso_string+'</EndTime>').firstChild
            output_dom.firstChild.appendChild( output_dom.createTextNode('\n'))
            output_dom.firstChild.appendChild(end_xml)

        if initialStepInfo.wall != None and finalStepInfo.wall != None:
            difference = int( float(finalStepInfo.wall) - float(initialStepInfo.wall) )
            hours = int( difference/3600 )
            minutes = int( (difference - hours*3600)/60 )
            runtime_xml = minidom.parseString('<RunTime>%d:%02d</RunTime>'%(hours,minutes)).firstChild
            output_dom.firstChild.appendChild( output_dom.createTextNode('\n'))
            output_dom.firstChild.appendChild(runtime_xml)
        sys.stdout.write( output_dom.toxml())
        sys.stdout.write("\n")

    else:
        outdict = {}
        outdict['dti'] = initialStepInfo.dt
        outdict['dtf'] = finalStepInfo.dt
        outdict['timei'] = initialStepInfo.time
        outdict['timef'] = finalStepInfo.time
        outdict['cyclei'] = initialStepInfo.cycle
        outdict['cyclef'] = finalStepInfo.cycle

        for line in last_ten_lines:
            print(line[:-1])
        print("\n\n")

        #print "%(dti)s %(dtf)s %(timei)s %(timef)s %(cyclei)s %(cyclef)s"%outdict

    return {'initialStepInfo':initialStepInfo,'finalStepInfo':finalStepInfo, 'all_steps':all_steps}




def parse_perf(initialStepInfo,finalStepInfo,fname='performance.out'):
    if len(glob.glob(fname)) == 0:
        print("NO PERFORAMNCE")
        return {'proc':-1,'coresec':-1,'cellup':-1, 'cspercu':-1}

    fptr = open(fname,'r')
    ncore_re = re.compile(r'# Starting performance log. MPI processes:\s*(\d*)')
    cycle_re = re.compile(r'Cycle_Number\s*(\d*)')
    taking_data = False
    total_time = 0
    total_up = 0
    for line in fptr:
        match = ncore_re.match(line)
        if match != None:
            ncore=int( match.group(1))
            continue
        match = cycle_re.match(line)
        if match is not None:
            cycle = int(match.group(1))
            if cycle >= initialStepInfo.cycle:
                taking_data=True
            if cycle > finalStepInfo.cycle:
                taking_data=False
            continue
        if line.startswith('Total'):
            bits = []
            for bit in line.split(" "):
                if len(bit):
                    #get rid of extra spaces
                    bits.append(bit)
            # mean time(1), stddev time(2), Tmin (3), Tmax(4), Nupdates (5), Ngrids (6), mean cell updates/s/processor(7)
            # I believe that max time is the one to take, as it sets the performance
            # for the whole level.
            total_time += float(bits[4])
            total_up += float(bits[5])
    return {'proc':ncore,'coresec':ncore*total_time,'cellup':total_up, 'cspercu':(ncore*total_time/total_up)}

all_stuff = None
plot_name = "tparse_multi"
base_name = []
print(args)
if len(args) > 1:
    for fname in args:
        base1 = fname.split('.')[0]
        if base1 not in base_name:
            base_name.append(base1)
            plot_name += "%s_"%base1
        all_stuff=scrub_log_file(fname,all_stuff)
else:
    fname = args[0]
    plot_name = 'tparse_'+fname
    all_stuff=scrub_log_file(fname,all_stuff)
        

if 0:
    print("No plots")
else:
    def add_dumps(plot_obj, full_list, in_or_out, cycle_or_time, y_value=-1,min_x=-1, print_number='some'):
        if print_number == 'no':
            print("HARD NO")
            pdb.set_trace()
        output_list = full_list[in_or_out]
        if len(output_list) == 0:
            return
        all_x = nar([ output[cycle_or_time] for output in output_list])
        c={'output':'b','input':'r'}[in_or_out]
        all_dumpnum = [output['number'] for output in output_list]
        all_x [ all_x < min_x ] = min_x
        plot_obj.scatter(all_x, np.zeros_like(all_x) + y_value,c=c)
        if print_number == 'some':
            plot_obj.text(all_x[0], y_value*1.1, "%s%04d"%(output_list[0]['dir'], all_dumpnum[0]),color=c)
            if len(all_x) > 1:
                plot_obj.text(all_x[-1], y_value*1.1, "%s%04d"%(output_list[-1]['dir'], all_dumpnum[-1]),color=c)
        elif print_number == 'all':
            for i in range(len(all_x)):
                plot_obj.text(all_x[i], y_value*1.1, "%s%04d"%(output_list[i]['dir'], all_dumpnum[i]),color=c)
    all_steps = all_stuff['all_steps']
    all_cycle  =nar(list(map(int, all_steps['cycle'])))
    all_dt     =nar(list(map(float,all_steps['dt'])))
    all_time   =nar(list(map(float,all_steps['time'])))
    plot_format = 'png'
    warning_check=False
    for scale in ['log','linear']:
        plt.clf()
        plt.plot(all_cycle, all_dt)
        add_dumps(plt,all_steps,'output', 'cycle',min_x=all_cycle.min(), y_value=all_dt.min())
        plt.xlabel('cycle'); plt.ylabel('dt')
        this_outname = plot_name+'cycle_dt.%s'%plot_format
        if scale == 'log':
            ylim = plt.ylim()
            if ylim[0] < 0:
                warning_check=True
            plt.yscale('log')
            add_dumps(plt,all_steps,'input', 'cycle',min_x=all_cycle.min(),y_value= all_dt.min()*1.5,print_number='all')
            plt.ylim( ylim)
            this_outname = plot_name+'cycle_dt_log.%s'%plot_format
        plt.savefig(this_outname)

    for scale in ['log','linear']:
        plt.clf()
        plt.plot(all_time, all_dt)
        add_dumps(plt,all_steps,'output', 'time',min_x=all_time.min(),y_value= all_dt.min(), print_number='some')
        plt.xlabel('time'); plt.ylabel('dt')
        this_outname = plot_name+'time_dt.%s'%plot_format
        if scale == 'log':
            ylim = plt.ylim()
            if ylim[0] < 0:
                warning_check=True
            add_dumps(plt,all_steps,'input', 'time',min_x=all_time.min(),y_value= all_dt.min()*2,print_number='all')
            plt.yscale('log')
            plt.ylim( ylim)
            this_outname = plot_name+'time_dt_log.%s'%plot_format
        plt.savefig(this_outname)
    if warning_check:
        if np.sum( np.isnan(all_dt) ) == 0 and all_dt.min() > 0:
            print("IGNORE THE SCALE WARNING.  dt is not nan and not negative.  This is the limits on the linear plot.")

#
#if 'all_steps' in log_out:
#    filename = 'butts'
#    all_steps = log_out['all_steps']
#    plt.clf()
#    plt.plot(all_steps['cycle'], all_steps['dt'])
#    plt.xlabel('n'); plt.ylabel('n')
#    plt.savefig(filename+'_cycle_dt.%s')
#   plt.yscale('log')
#   plt.savefig(filename+'_cycle_dt_log.%s')
#   plt.clf()
#   plt.plot(all_steps['time'], all_steps['dt'])
#   plt.xlabel('time'); plt.ylabel('dt')
#   plt.savefig(filename+'_time_dt.%s')
#   plt.yscale('log')
#   plt.savefig(filename+'_time_dt_log.%s')
#dumb_Initial=stepInfo
#dumb_Final = stepInfo
#dumb_Initial.cycle = 0
#dumb_Final.cycle = 5
#perf_out=parse_perf(dumb_Initial,dumb_Final)
perf_out = parse_perf(all_stuff['initialStepInfo'],all_stuff['finalStepInfo'])
outdict = {}
outdict['dti'] =   all_stuff['initialStepInfo'].dt
outdict['dtf'] =   all_stuff['finalStepInfo'].dt
outdict['timei'] = all_stuff['initialStepInfo'].time
outdict['timef'] = all_stuff['finalStepInfo'].time
outdict['cyclei'] =all_stuff['initialStepInfo'].cycle
outdict['cyclef'] =all_stuff['finalStepInfo'].cycle
outdict.update(perf_out)

print("%(dti)s %(dtf)s %(timei)s %(timef)s %(cyclei)s %(cyclef)s %(cellup)0.3e %(coresec)0.2e %(cspercu)0.3e"%outdict)
print("then also", perf_out)


#end
