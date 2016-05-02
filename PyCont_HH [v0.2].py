""" EXAMPLE: Hodgkin-Huxley
    Kryk, April 2016
"""

from PyDSTool import *

vfn_str = '(Iapp-ionic(v,m,h,n))/C'
mfn_str = 'ma(v)*(1-m)-mb(v)*m'
nfn_str = 'na(v)*(1-n)-nb(v)*n'
hfn_str = 'ha(v)*(1-h)-hb(v)*h'
aux_str = 'm*m*m*h'

#Parameters from Guckenhemer 1993(Bifurcation of the Hodgkin-Huxley...)
pars_args = {
        'Iapp': -115,
        'C': 1.0,

        'vK': 12.0,
        'gK': 12.0,

        'vNa': -115.0,
        'gNa': 120.0,

        'vL': 10.95,
        'gL': 0.3,}

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
DSargs.fnspecs = auxdict
DSargs.ics = ic_args

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

#----------------------------------
Iapp = DSargs.pars['Iapp']
gK = DSargs.pars['gK']
gNa = DSargs.pars['gNa']
gL = DSargs.pars['gL']

max_Iapp = 120

default_gK = 12
step_gK = -10
max_gK = 42

default_gNa = 120
step_gNa = -40
max_gNa = 200

default_gL = 0.3
step_gL = 1.
max_gL = 3
#---------------------------------

while Iapp < max_Iapp:

	start = clock()
	PyCont['EQ1'].forward()
	PyCont['EQ1'].backward()

	# Plot
	PyCont.display(('Iapp','v'),stability=True)
	#PyCont['LC1'].display(('Iapp','v_min'),stability=True)
	PyCont['EQ1'].display(coords=('Iapp', 'v'), stability=True)

	plt.xlim([-400, 400])
	PyCont.plot.fig1.axes1.axes.set_title('Bifurcation Diagram')
	PyCont.plot.fig1.toggleAll('off', bytype=['P', 'MX'])

	#PyCont.plot.setLegends('_nolegend_', bytype=['P', 'MX'])
	#PyCont.plot.fig2.refresh()
	#plt.legend()

	
	plt.savefig('C:\Users\Smash\Pictures\Test\(Iapp-%d,gK-%d,gNa-%d,gL-%d).png' % (Iapp,gK,gNa,gL))
	plt.close()

	if gK < max_gK and gNa == default_gNa and gL == default_gL:
		gK = gK + step_gK
	elif gK == max_gK and gNa < max_gNa:
		gK = default_gK
		gNa = gNa + step_gNa


	elif gK == default_gK and gNa > default_gNa and gNa < max_gNa:
		gNa = gNa + step_gNa
	elif gNa == max_gNa and gL < max_gL:
		gNa = default_gNa
		gL = gL + step_gL


	elif gNa == default_gNa and gL > default_gL and gL < max_gL:
		gL = gL + step_gL	


	elif gL == max_gL and Iapp == -210:
		Iapp = -115
		gL = default_gL
	elif gL == max_gL and Iapp == -115:
		Iapp = -50
		gL = default_gL
	elif gL == max_gL and Iapp == -50:
		Iapp = 110
		gL = default_gL
	elif gL == max_gL and Iapp == 110:
		Iapp = max_Iapp