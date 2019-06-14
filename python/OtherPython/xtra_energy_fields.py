from yt import ValidateParameter, add_field
import numpy as np
field_args={}
def dave_add_field(yt_object, field_name=None, argument_dict=field_args):
    if field_name is not None:
        field_list = [field_name]
    else:

        field_list = ['momentum_x' ,'momentum_y' ,'momentum_z' ,\
                      'momentum_magnitude' ,'mean_square_velocity' ,\
                      'eng_x' ,'eng_y' ,'eng_z' ,'rel_kinetic_energy' ,\
                      'grav_pot' ,'gas_work']
    for field in field_list:
        yt_object.add_field(field,sampling_type='cell',**argument_dict[field])

momentum_units = 'g/(cm**2*s)'
def _momentum_x(field,data):
    return data['density']*data['velocity_x']
field_args['momentum_x']={'function':_momentum_x,'units':momentum_units}
def _momentum_y(field,data):
    return data['density']*data['velocity_y']
field_args['momentum_y']={'function':_momentum_y,'units':momentum_units}
def _momentum_z(field,data):
    return data['density']*data['velocity_z']
field_args['momentum_z']={'function':_momentum_z,'units':momentum_units}
def _momentum_magnitude(field,data):
    out = data['momentum_x']**2+data['momentum_y']**2 + data['momentum_z']**2
    return np.sqrt(out)
field_args['momentum_magnitude']={'function':_momentum_magnitude,'units':momentum_units}

def _mean_square_velocity(field,data):
    out = data['velocity_x']**2+data['velocity_y']**2 + data['velocity_z']**2
    return out
field_args['mean_square_velocity']={'function':_mean_square_velocity,'units':"cm**2/s**2"}

eng_units = 'g/(cm*s**2)'
def _eng_x(field,data):
    return 0.5*data['density']*data['velocity_x']*data['velocity_x']
field_args['eng_x']={'function':_eng_x,'units':eng_units}
def _eng_y(field,data):
    return 0.5*data['density']*data['velocity_y']*data['velocity_y']
field_args['eng_y']={'function':_eng_y,'units':eng_units}
def _eng_z(field,data):
    return 0.5*data['density']*data['velocity_z']*data['velocity_z']
field_args['eng_z']={'function':_eng_z,'units':eng_units}

def _rel_kinetic_energy(field,data):
    if data.has_field_parameter('bulk_velocity'):
        vbar = data.get_field_parameter('bulk_velocity')
    else:
        vbar = data.ds.arr([0]*3,'code_velocity')


    vx = data['velocity_x']-vbar[0]
    vy = data['velocity_y']-vbar[1]
    vz = data['velocity_z']-vbar[2]
    return 0.5*data['density']*(vx*vx+vy*vy+vz*vz)
field_args['rel_kinetic_energy']={'function':_rel_kinetic_energy,'units':eng_units,'validators':[ValidateParameter('bulk_velocity')]}

def _grav_pot_grad(field,data):
    gx = grad(data,'PotentialField',0)
    gy = grad(data,'PotentialField',1)
    gz = grad(data,'PotentialField',2)

def _grav_pot(field,data):
    try:
        output = 0.5*data['density']*data['PotentialField']
    except:
        output = data['density']*0

    return output
field_args['grav_pot']={'function':_grav_pot,'units':eng_units}

def _gas_work(field,data):
    return (data['density'].in_units('code_density').v*np.log(data['density'].in_units('code_density').v))*data.ds.quan(1,eng_units)
field_args['gas_work']={'function':_gas_work,'units':eng_units}


