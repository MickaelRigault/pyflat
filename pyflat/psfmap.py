import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.stats import multivariate_normal


class PSFMap(object):

	def __init__(self,width,height):
		self.width=width
		self.height=height
		self.map=np.zeros((self.width+202,self.height+202))+10**-200

	def get_psf_map(self,x,y,f,seeing):
		""" get a array (width x height) with simulated stars located at x, y
			properties of the stars are the sep flux and the seeing of the image
	
		"""  
		#input correction allow user to input floats instead of arrays
		if not hasattr(seeing,'__iter__'):
			seeing=[seeing]

		if not hasattr(x,'__iter__'):
			x=[x]

		if not hasattr(y,'__iter__'):
			y=[y]

		if not hasattr(f,'__iter__'):
			f=[f]


		if len(x)!=len(y):
			sys.exit('x and y must be of same length')
		else:
			if len(seeing)!=len(x):
				if len(seeing)!=1:
					sys.exit('seeing and x,y must be of same length, or of length 1')
				else:
					seeing=np.ones(len(x))*seeing

				if np.max(x)>self.width:
					sys.exit('check image orientation')
				if np.max(y)>self.height:
					sys.exit('check image orientation')
							
				for i in range(0,len(x)):
					if i%100==0:
						print i, '/', len(x)
					self.map[int(x[i]+1):int(x[i]+201),int(y[i]+1):int(y[i]+201)]+=self.get_psf(x[i],y[i],seeing[i],200,'2gauss')*f[i]
	
		return self.map[101:-101,101:-101]
	
	def get_psf(self,a,b,c,d,param):
		if param=='2gauss':
			return self._psf_2gaussians(a,b,c,d) 

	def show(self):
		plt.imshow(np.log10(self.map[101:-101,101:-101]),interpolation='nearest',vmin=-15)
		plt.show()

	def _psf_2gaussians(self,xp,yp,seeing,w):
		xg,yg=np.mgrid[0:w:1,0:w:1]
		pos=np.empty(xg.shape+(2,))
		pos[:, :, 0]=xg; pos[:, :, 1]=yg
		gauss1=multivariate_normal.pdf(pos,mean=[w/2.+(xp%1)-0.01386,w/2.+(yp%1)-0.03465],cov=[[8/9.,.0],[.0,8/9.]])
		gauss2=multivariate_normal.pdf(pos,mean=[w/2.+(xp%1)-0.24257,w/2.+(yp%1)-0.62871], cov=[[5.,.0],[.0,5.]])
		return (gauss1*0.9+gauss2*0.1)*0.877857887





