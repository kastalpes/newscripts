#from go import *
from GL import *
import matplotlib.colors as colors
def plotter2(arrays,fname,norm=None,labels=['Q','U','E','B'],axis_labels=None,npx=2,
             share=True,zmin=None,**args):
    print('WA',labels)
    np = len(arrays)
    npy = np//npx
    if np%npx: npy+=1
    fig, axes = plt.subplots(npx,npy,sharex=share,sharey=share, figsize=(12,12))
    print("NPX %d NPY %d"%(npx,npy))
    max_val = max([arr.max() for arr in arrays])
    min_all = min([arr.min() for arr in arrays])
    if zmin is not None:
        min_all = zmin
    if norm is 'positive':
        min_val = min([ arr[arr>0].min() for arr in arrays])
        args['norm']=colors.LogNorm(vmin=min_val,vmax=max_val)
        print("WTF",min_val)
    elif norm is 'symlog' :
        min_val = min_all
        args['norm']=colors.SymLogNorm(linthresh=1e-10,vmin=min_val,vmax=max_val)
    elif norm is 'ind':
        args['norm']=None
    else:
        min_val = min_all
        args['norm']=colors.Normalize(vmin=min_val,vmax=max_val)

    args['interpolation']='nearest'
    args['origin']='lower'
    oot=[]
    for n,arr in enumerate(arrays):
        ipy = n//npx
        ipx = n%npx
        if npx == 1:
            this_ax = axes
        elif npy == 1:
            this_ax = axes[ipx]
        else:
            this_ax = axes[ipx][ipy]
        oot.append(this_ax.imshow(arr,**args))
        if norm is 'ind': cb=fig.colorbar(oot[-1],ax=this_ax)
        this_ax.set_title(labels[n])
#   oot.append(axes[1][0].imshow(U,**args))
#   if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][0])
#   axes[1][0].set_title(labels[1])
#   oot.append(axes[0][1].imshow(E,**args))
#   if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[0][1])

#   axes[0][1].set_title(labels[2])
#   oot.append(axes[1][1].imshow(B,**args))
#   if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][1])
#   axes[1][1].set_title(labels[3])
    if norm is not 'ind':
        cb=fig.colorbar(oot[-1],ax=axes)
        if norm is 'positive':
            cb.cmap.set_under('w')
    if axis_labels:
        for n in range(np):
            axes[n//2][n%2].set_xlabel(axis_labels[0])
            axes[n//2][n%2].set_ylabel(axis_labels[1])
    fig.savefig(fname)
    print(fname)
    plt.close(fig)

def plotter(Q,U,E,B,fname,norm=None,labels=['Q','U','E','B'],axis_labels=None,**args):
    fig, axes = plt.subplots(2,2,sharex=True,sharey=True)
    max_val = max([Q.max(),U.max(),E.max(),B.max()])
    if norm is 'positive':
        min_val = 1e-15
        args['norm']=colors.LogNorm(vmin=min_val,vmax=max_val)
    elif norm is 'symlog' :
        min_val = min([Q.min(),U.min(),E.min(),B.min()])
        args['norm']=colors.SymLogNorm(linthresh=1e-10,vmin=min_val,vmax=max_val)
    elif norm is 'ind':
        args['norm']=None
    else:
        min_val = min([Q.min(),U.min(),E.min(),B.min()])
        args['norm']=colors.Normalize(vmin=min_val,vmax=max_val)

    args['interpolation']='nearest'
    args['origin']='lower'
    oot=[]
    oot.append(axes[0][0].imshow(Q,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[0][0])
    axes[0][0].set_title(labels[0])
    oot.append(axes[1][0].imshow(U,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][0])
    axes[1][0].set_title(labels[1])
    oot.append(axes[0][1].imshow(E,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[0][1])

    axes[0][1].set_title(labels[2])
    oot.append(axes[1][1].imshow(B,**args))
    if norm is 'ind': cb=fig.colorbar(oot[-1],ax=axes[1][1])
    axes[1][1].set_title(labels[3])
    if norm is not 'ind':
        cb=fig.colorbar(oot[-1],ax=axes)
        if norm is 'positive':
            cb.cmap.set_under('w')
    if axis_labels:
        for n in range(4):
            axes[n//2][n%2].set_xlabel(axis_labels[0])
            axes[n//2][n%2].set_ylabel(axis_labels[1])
    fig.savefig(fname)
    plt.close(fig)
