'''
Author: Charalampos Mavidis 
E-Mail: mavidis@iesl.forth.gr
TransferMatrix Method

			 x
			 ^
			 |
eps1,mu1	 |   eps2, mu2
			 |
_____________|__________________>z
			 |
			 |
             |
             sigma


'''

'''
Convension

(hin) = (M)(hout)
'''


import numpy as np
import math
from scipy import constants
import copy

eps0 = constants.epsilon_0
mu0 = constants.mu_0
c = constants.c
''' P-polarization

E = (Ex,0,Ez) ; H = (0,H,0)

'''
def M_p(freq, k1z, k2z, eps1, eps2):
	
	eta = (eps1/eps2)*(k2z/k1z)

	M11 = 1+eta + 0j
	M12 = 1-eta + 0j
	M21 = 1-eta + 0j
	M22 = 1+eta + 0j
	M = 0.5*np.matrix([[M11,M12],
					   [M21,M22]])
	return M

''' S-Polarization

E = (0,E,0) ; H = (Hx,0,Hz)

 '''
def M_s(freq, k1z, k2z, mu1, mu2):
	xi  = mu1*mu0*freq*sigma/k1z
	eta = mu1*k2z/(mu2*k1z)

	zeta = k2z*sigma_m/(freq*mu2*mu0)
	#print eta,xi

	M11 = 1.0+eta + 0j
	M12 = 1.0+eta + 0j
	M21 = 1.0-eta + 0j
	M22 = 1.0-eta + 0j
	M = 0.5*np.matrix([[M11,M12],
					   [M21,M22]])
	return M

''' Propagation Matrix '''
def Propagate(k,d):
	a = np.exp(-1j*k*d) + 0j
	b = np.exp(+1j*k*d) + 0j
	P = np.matrix([[a,0],
				   [0,b]])
	return P


''' Construct Transfer Matrix for a given list of layers ''' 
class TransferMatrix(object): 

	def __init__(self,structure):
		"""Construct TransferMatrix
		
		Arguments:
			structure {[list]} -- List of layers 
		"""
		self.Ms = 1.0+0j
		self.Mp = 1.0+0j

		self.ts = 0+0j
		self.rs = 0+0j
		self.rp = 0+0j
		self.tp = 0+0j

		self._runSetUp = False

		# Check if structure is a list
		if not type(structure) is list: 
			assert "A list of layers is required."
			return -1 

		self.structure = structure

		# Construct the Transfer Matrices
		self.setup()

	def setup(self): 
		# Need an initial layer: layer0
		layer0 = self.structure[0]
		for layer1 in self.structure[1:]: 
			if isinstance(layer1, Layer): 
				self.Ms *= M_s(layer1.freq, layer0.kz, layer1.kz, layer0.kz, layer1.kz, layer1.sigma,sigma_m=0)
				self.Ms *= Propagate(layer1, layer1.thickness)

				self.Mp *= M_s(layer1.freq, layer0.kz, layer1.kz, layer0.eps, layer1.eps, layer1.sigma,sigma_m=0)
				self.Mp *= Propagate(layer1.kz, layer1.thickness)

				# Next layer
				layer0 = layer1
		self._runSetUp = True
 
	def calculate(self): 
		'''
		Calculate transmission and reflection coefficients
		''' 
		if self._runSetUp:
			self.ts = 1./self.Ms[0,0] + 0j
			self.rs = self.Ms[1,0]/Ms[0,0] + 0j
			self.rp = 1./self.Mp[0,0] + 0j
			self.rp = self.Mp[1,0]/Mp[0,0] + 0j





class Layer(object):

	def __init__(self,eps,mu,sigma,d,kz,freq,name=""):
		self.epsilon = eps
		self.n = np.sqrt(eps*mu)
		self.sigma = sigma
		self.mu = mu
		self.name = name
		self.thickness = d
		self.d = d # just for backward compatibility

		self.kz = kz



		if abs(np.real(self.n)) > 0:
			self.wavelength = constants.c/freq/np.real(self.n)
		else:
			self.wavelength = constants.c/freq/100.0

		# Number of cells
		self.Dz = self.thickness/200.0
		self.Nm = math.ceil(self.thickness/self.Dz)

		self.position = np.linspace(0,d,self.Nm)


		# Field at the interface between this layer and layer i+1
		self.field = np.matrix([[0],[0]],dtype=np.complex)
		# Propagate the field
		self.propagated_field = np.matrix([[0],[0]],dtype=np.complex)

		self.field_array = np.zeros(self.Nm,dtype=np.complex)

	def get_field(self,field_type="complex"):
		if field_type=="real":
			return np.real(self.field_array[::-1])
		elif field_type=="imag":
			return np.imag(self.field_array[::-1])
		elif field_type=="abs":
			return np.abs(self.field_array[::-1])
		elif field_type=="phase":
			return np.angle(self.field_array[::-1])
		else:
			return self.field_array[::-1] # complex


class SeminfiniteLayer(Layer):

	def __init__(self,eps,kz,freq):
		super(SeminfiniteLayer, self).__init__(eps,1,0,0,kz,freq,name="")

 
''' Field calculation'''
class Field(object):

	def __init__(self,t,r,freq,polarization):
		'''
		Field calculation
		Parameters needed: t, r, structure
		Need M11*t,M21*t
		'''
		self._layers = list()
		self._freq = freq

		self.total_thickness = 0.0

		self.transmitted_field = np.matrix([[t],[0]],dtype=np.complex)
		self.t = t
		self.r = r
		self.reflected = np.matrix([[0],[r]],dtype=np.complex)
		self.polarization = polarization

		self.number_of_cells = 0
		self.substrate = None
		self.superstrate_inc  = None
		self.superstrate_r  = None

		self.structure = []

	def __call__(self): 
		for lay in self._layers:
			print lay.name, " | "


	def add_substrate(self,lay):
		self.substrate = lay

	def add_superstrate(self,lay):
		self.superstrate = lay # total field
		self.superstrate_r = copy.deepcopy(lay) # reflected field
		self.superstrate_inc = copy.deepcopy(lay) # incident field


	def add_layer(self,lay):
		'''
		add_layer:
		@param lay: Layer object

		Note: The layers should be added in reverse order.
		'''
		self._layers.append(lay)
		self.total_thickness += lay.d

	def setup(self):
		'''
		setup()
		Prepare position arrays
		'''
		z_start = 0.0
		z_end = self._layers[0].d
		for lay in self._layers:
			lay.position = np.linspace(z_start,z_end,lay.Nm)
			z_start = z_end
			z_end = z_start + lay.d

			self.structure.append(np.real(lay.epsilon)*np.ones(lay.Nm).flatten().tolist())
		self.structure = np.array(self.structure).flatten()


	
	def calculate(self):
		'''
		At least one layer has to be in self._layers
		'''
		if len(self._layers)==0:
			assert "No layers."
			return -1

		self.setup()
		# Initialize

		if self.polarization in ['s','p']:

			# Propagate in substate layer
			self.substrate.field_array[self.substrate.Nm-1] = self.t
			fl = self.t
			for cell_index in np.arange(self.substrate.Nm-2,-1,-1):
				self.substrate.field_array[cell_index] = fl*np.exp(1j*self.substrate.kz*self.substrate.Dz)
				fl = self.substrate.field_array[cell_index]

			self.substrate.propagated_field = self.transmitted_field
			forward_layer = self.substrate

			for lay in self._layers:
					# M_s(freq, k1z, k2z, mu1, mu2, sigma)
					if self.polarization == 's': 
						tr_matrix = M_s(self._freq, lay.kz, forward_layer.kz, lay.mu, forward_layer.mu, forward_layer.sigma)
					elif self.polarization == 'p':
						tr_matrix = M_p(self._freq, lay.kz, forward_layer.kz, lay.epsilon, forward_layer.epsilon, forward_layer.sigma)
					lay.field = tr_matrix*forward_layer.propagated_field

					fll =  Propagate(lay.kz,lay.d)*lay.field

					print np.real(np.sum(fll)),np.real(np.sum(lay.field)),self.t

					# Change previous layer
					forward_layer = lay

					# Fill in the field array
					lay.field_array[0] = np.sum(lay.field)
					field_i = lay.field
					#print lay.field, np.sum(lay.field)
					for cell_index in np.arange(1,lay.Nm):
						field_i = Propagate(lay.kz,lay.Dz)*field_i
						#print field_i
						lay.field_array[cell_index] = np.sum(field_i)
					# Propagate back
					lay.propagated_field = field_i
					#print field_i.shape,field_i
					self.forward_layer = lay

			if self.polarization == 's': 
				tr_matrix = M_s(self._freq, self.superstrate.kz, forward_layer.kz, self.superstrate.mu, forward_layer.mu, lay.sigma)
			elif self.polarization == 'p':
				tr_matrix = M_p(self._freq, self.superstrate.kz, forward_layer.kz, self.superstrate.epsilon, forward_layer.epsilon, lay.sigma)
			field_i = tr_matrix*forward_layer.propagated_field
			#print field_i
			# Propagate in the superstrate layer

			#print self.r, np.real(np.sum(self.superstrate.field))
			#self.superstrate_inc.field_array[self.superstrate_inc.Nm-1] = 1
			flr = self.r #lay.field_array[cell_index]-1.0#self.r
			flinc = 1.0

			self.superstrate.field_array[self.superstrate.Nm-1] = np.sum(field_i) # total field
			self.superstrate_r.field_array[self.superstrate_r.Nm-1] = self.r # reflected field
			self.superstrate_inc.field_array[self.superstrate_inc.Nm-1] = 1.0 # total field
			for cell_index in np.arange(self.superstrate_r.Nm-2,-1,-1):
				self.superstrate_r.field_array[cell_index] = flr*np.exp(+1j*self.superstrate_r.kz*self.superstrate_r.Dz)
				flr = self.superstrate_r.field_array[cell_index]

				self.superstrate_inc.field_array[cell_index] = flinc*np.exp(-1j*self.superstrate_inc.kz*self.superstrate_r.Dz)
				flinc = self.superstrate_inc.field_array[cell_index]

				field_i = Propagate(self.superstrate.kz,self.superstrate.Dz)*field_i
				
				self.superstrate.field_array[cell_index] = np.sum(field_i)
				#flinc = self.superstrate_inc.field_array[cell_index]
