""" EXAMPLE: Hodgkin-Huxley
    Bifircation Diagram (Iapp, V)
    Kryk, April 2016
"""

from PyDSTool import *
plt.ion()

vfn_str = '(Iapp-ionic(v,m,h,n))/C'
mfn_str = 'ma(v)*(1-m)-mb(v)*m'
nfn_str = 'na(v)*(1-n)-nb(v)*n'
hfn_str = 'ha(v)*(1-h)-hb(v)*h'
aux_str = 'm*m*m*h'

#Parameters from Guckenhemer 1993(Bifurcation of the Hodgkin-Huxley...)
pars_args = {
        'Iapp': 10.0,
        'C': 1.0,

        'vK': -100.0,
        'gK': 80.0,

        'vNa': 100.0,
        'gNa': 50.0,

        'vL': -67.0,
        'gL': 0.1, }

ic_args = {
        'v':-70.0, 
        'm': 0, 
        'h': 1, 
        'n': 0 }

auxdict = {
        'ionic': (['v', 'm', 'h', 'n'], 
    'gNa*m*m*m*h*(v-vNa) + gK*n*n*n*n*(v-vK) + gL*(v-vL)'),

        'ma': (['v'], '0.32*(v+54)/(1-exp(-(v+54)/4))'),
        'mb': (['v'], '0.28*(v+27)/(exp((v+27)/5)-1)'),
        'ha': (['v'], '.128*exp(-(50+v)/18)'),
        'hb': (['v'], '4/(1+exp(-(v+27)/5))'),
        'na': (['v'], '.032*(v+52)/(1-exp(-(v+52)/5))'),
        'nb': (['v'], '.5*exp(-(57+v)/40)') }

DSargs = args(name='HH')
DSargs.pars = pars_args
DSargs.varspecs = {
                    'v': vfn_str, 'm': mfn_str,
                    'h': hfn_str, 'n': nfn_str }
                    #,
                    #'v_bd0': 'getbound("v",0)',
                    #'v_bd1': 'getbound("v",1)' }
DSargs.fnspecs = auxdict
DSargs.ics = ic_args
DSargs.xdomain = {'v': [-130, 70], 'm': [0,1], 'h': [0,1], 'n': [0,1]}

testDS = Generator.Vode_ODEsystem(DSargs)

# Set up continuation class
PyCont = ContClass(testDS)

PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['Iapp']
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 350
PCargs.MaxStepSize = 1.
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PyCont.newCurve(PCargs)

print('Computing curve...')
start = clock()
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()
print('done in %.3f seconds!' % (clock()-start))

# Plot
PyCont.display(('Iapp','v'),stability=True)
#PyCont['LC1'].display(('Iapp','v_min'),stability=True)
PyCont['EQ1'].display(coords=('Iapp', 'v'), stability=True)

plt.xlim([0, 600])
PyCont.plot.fig1.axes1.axes.set_title('Bifurcation Diagram')
PyCont.plot.fig1.toggleAll('off', bytype=['P', 'MX'])

#PyCont.plot.setLegends('_nolegend_', bytype=['P', 'MX'])
#PyCont.plot.fig2.refresh()
#plt.legend()
show()