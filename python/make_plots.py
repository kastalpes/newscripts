
from GL import *
import plot_quantities as pq
reload(pq)

out_prefix = 'ze02'
directory = '/scratch/00369/tg456484/Paper49d_moresims/ze01_M10_MA1_1024_quan'

out_prefix='za02'
directory='/scratch/00369/tg456484/Paper49d_moresims/za02_M10_MA1_64'

out_prefix='zb02'
directory='/scratch/00369/tg456484/Paper49d_moresims/zb02_M10_MA1_128'

out_prefix='zc02'
directory='/scratch/00369/tg456484/Paper49d_moresims/zc02_M10_MA1_256'

out_prefix='zd02'
directory='/scratch/00369/tg456484/Paper49d_moresims/zd02_M10_MA1_512'
out_prefix='za03'
directory='/scratch/00369/tg456484/Paper49d_moresims/za03_goddamnit_64/'
out_prefix='zb03'
directory='/scratch/00369/tg456484/Paper49d_moresims/zb03_goddamnit_128'

pq.v2_b2_plot(out_prefix,directory)
