from flask import Flask, request, render_template ,  Response, jsonify
import matplotlib.pyplot as plt
import math
import io
import base64

import math     as mt
import numpy	as np

import datetime as dt

class   Outrigger:
    
    def __init__(self, data):
        self.Htot = float ( data['Htot'])
        self.D = float ( data['D'])
        self.B = float ( data['B'])
        self.Hstory = float ( data['Hstory'])
        self.Load = float ( data['Load'])
        
        self.EGirder =float (  data['EGirder'])
        self.ECore =float (  data['ECore'])
        self.Ecol = float ( data['Ecol'])

        self.IGirder = float ( data['IGirder'])
        self.ICore = float ( data['ICore'])
        self.ACol = float ( data['ACol'])

        self.Nout = int ( data['Nout'])
        
        #x = np.zeros (self.Nout)
        #for ii in range (0 , self.Nout,1):
         #   self.x[ii] = float(data['x'][ii])

        xx = data ['x']
        
        self.x = [float(num) for num in xx.split(",")]
        print (self.Nout , self.x)


    def Out_Analysis   (self):
        n       = self.Nout
        a       = self.D/2-self.B  
        EIcore  = self.ICore * self.ECore
        EIGirder= self.IGirder * (1+a/self.B)**3 * self.EGirder    
        EACol   = self.ACol * self.Ecol
        
        x1 = np.zeros (n)
        A = np.zeros ((n,n))
        eta = np.zeros (n)


    #           Constants 
        beta  = EIcore / EIGirder  *self.D/self.Htot;
        alpha = EIcore /( EACol * (self.D*self.D/2));

        omega = beta  /(12* (1+alpha));
        s = 1/EIcore  + 2/ ( self.D * self.D * EACol);

        fact = self.Load    *self.Htot*self.Htot / (6*EIcore * s) 
        
        for ii in range (0 , n,1):
            eta [ii]= self.x[ii]/ self.Htot
          
        for ii in range (0 , n,1):
            A[ii,ii] = omega + (1-eta[ii]);
            x1[ii] = (1-eta[ii]**3) *fact;
        
            for iii in range  (ii+1,n,1):							
                A[ii,iii] =  (1-eta[iii]);							3	 				 
                A[iii,ii] = A[ii,iii]
    #                               solve the Equations    
      
        M = np.linalg.solve(A, x1)

        Nst = int( self.Htot/ self.Hstory)

        self.MwC = np.zeros (Nst+1)
        self.TcC = np.zeros (Nst+1)
        self.DwC = np.zeros (Nst+1)
       
        for  iii in range (0,Nst+1,1 ):	
            xh = iii*self.Hstory
            z = xh/self.Htot
            self.MwC[iii] = -self.Load*(xh*xh)/2             
            self.DwC[iii] = self.Load * self.Htot**4 / (24 * EIcore )  * (3-4*z+ z**3)
          #  self.Dw[iii] = self.DwC[iii]

        self.Mw = np.zeros (Nst+1)
        self.Tc = np.zeros (Nst+1)
        self.Dw = np.zeros (Nst+1)

        FlO =np.zeros (n+1,dtype=int)
        # self.Res    =   np.zeros((Nst+1,7),dtype=float)

        for  iii in range (0,n ):	
            FlO [iii] = int (self.x[iii]/self.Hstory)

        FlO[n]  =  Nst

        
        Macc = 0
        for iiii in range (0 , FlO[1]+1):                
                self.Mw[iiii] = self.MwC[iiii] 

        for  iii in range (0,n,1 ): 
            Macc = Macc + M[iii]                
            for iiii in range (FlO [iii] , FlO[iii+1]+1):                
                self.Mw[iiii] = self.MwC[iiii] + Macc 
                self.Tc[iiii] = self.TcC[iiii] + Macc / (self.D)


        for  iii in range (0,n,1 ):                         
            for iiii in range (0 , Nst+1 ):
                xh  =   iiii*self.Hstory
                xb  =   self.Htot-  xh
                z = xh/self.Htot
                if  (xh > self.x[iii] ):            # down
                    self.Dw[iiii]  = self.Dw[iiii] + M[iii] *xb**2 /(2*EIcore )
                   # print (xh,xb,self.Dw[iiii])
                else:                               # Up
                    self.Dw[iiii]  = self.Dw[iiii] +  M[iii] *(self.Htot -  self.x[iii])**2 /(2*EIcore )+ M[iii] *(self.Htot -  self.x[iii]) *(xb -(self.Htot-self.x[iii] ) )/EIcore 
                   # print (xh,xb,self.Dw[iiii])

        for iiii in range (0 , Nst+1 ):
            self.Dw[iiii] = self.DwC[iiii] - self.Dw[iiii]


        # output method 1              .
        #                                          RESULTS

        # self.Res =  "<link href={{ url_for('static', filename='stylesTHS.css') }} type='text/css' rel='stylesheet' />"
        
        self.Res =[]
        rec =  {'hei' :'Height' , 'axial': 'Wall Axial Force' ,'mom' :'Wall Moment' ,'defl' : 'Deflection'}
        self.Res.append (rec)


        for ix in range (0,Nst+1):
            hc = str( self.Htot-  ix * self.Hstory)
            rec =  {'hei' : hc ,  'axial':  "{:.2f}".format(self.Tc[ix] ),'mom' :"{:.2f}".format( self.Mw[ix]) ,'defl' : "{:.2f}".format( self.Dw[ix]*1000)}
            self.Res.append (rec)


            # self.Res  [hc ] = [ self.Tc[ix] , self.Mw[ix] , self.Dw[ix]]

            # self.Res  = {'hei' : hc}
            # self.Res  = {'axial': self.Tc[ix] }
            # self.Res  = {'mom' : self.Mw[ix] }
            # self.Res  = {'defl' : self.Dw[ix]}            
        
        rec =  {'hei' : "Cantilever" , 'axial': 0.00 ,'mom' :"{:.2f}".format( self.MwC[Nst] ) ,'defl' : "{:.2f}".format(self.DwC[0] * 1000)}
        self.Res.append (rec)
        # self.Res ["Cantilever" ] = [  0    ,self.MwC[Nst]  , self.DwC[0] * 1000  


        # self.Res ={'hei' : "Cantilever"}
        # self.Res = {'axial' :  0 }
        # self.Res = {'mom' : self.MwC[Nst] }
        # self.Res = {'defl' : self.DwC[0] * 1000  }

        return self.Res

app = Flask(__name__)

@app.route  ('/')

def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])


def analyze():
    My_Data = request.form.to_dict()
    outrigger = Outrigger(My_Data)
    result = outrigger.Out_Analysis()           
    plot1 = Draw_Geometry ( outrigger )       
    plot2 = Draw_Axial (outrigger)        
    plot3 = Draw_Moment ( outrigger )       
    plot4 = Draw_Deflection ( outrigger )       
    return render_template("index.html", plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4, Data=result)
    



def Draw_Geometry(outrigger):
        
    fig, ax = plt.subplots(figsize=(4, 5))   
        
    ax.plot([5 , 5  ], [ 0 ,outrigger.Htot ],  color ='black')  
    ax.plot([5+outrigger.D ,5+outrigger.D  ], [ 0 ,outrigger.Htot ], color ='black')  


    ax.plot([ 5+outrigger.B ,5+outrigger.B , 5+outrigger.D-outrigger.B , 5+outrigger.D-outrigger.B ], [ 0 ,outrigger.Htot   , outrigger.Htot  , 0] , color ='black')                  
    
    ax.plot([0 , 10+outrigger.D   ], [ 0  , 0], color ='black')  

    for ii in range (0, outrigger.Nout):

        ax.plot([5, 5+outrigger.B ], [ outrigger.Htot   - outrigger.x[ii] -2  ,  outrigger.Htot-  outrigger.x[ii] -2], color ='red')  
        ax.plot([5+outrigger.D - outrigger.B , 5+outrigger.D ], [ outrigger.Htot   - outrigger.x[ii] -2  ,  outrigger.Htot-  outrigger.x[ii] -2], color ='red')  
        ax.plot([5, 5+outrigger.B ], [ outrigger.Htot   - outrigger.x[ii] +2  ,  outrigger.Htot-  outrigger.x[ii] +2], color ='red')  
        ax.plot([5+outrigger.D - outrigger.B , 5+outrigger.D ], [ outrigger.Htot   - outrigger.x[ii] +2  ,  outrigger.Htot-  outrigger.x[ii] +2], color ='red')  

        # ax.plot([5+outrigger.B , 5+outrigger.D-outrigger.B ], [ outrigger.Htot   - outrigger.x[ii] -2  ,  outrigger.Htot-  outrigger.x[ii] -2])  
        # ax.plot([5+outrigger.B , 5+outrigger.D-outrigger.B ], [outrigger.Htot   - outrigger.x[ii] +2  ,  outrigger.Htot-  outrigger.x[ii] +2])  
        # ax.plot([0 , 10+outrigger.D   ], [ 0  , 0])  

    ax.set_title('Geometry') 
    ax.set_frame_on(True)

    ax.set_aspect(0.7)
    ax.set_axis_off ()
    

    # Save the plot to a BytesIO buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)

    # Serve the image as a response
    # return Response(buf, mimetype='image/png')
    return base64.b64encode(buf.read()).decode("utf-8")    

def Draw_Deflection( outrigger):
    yVal = np.arange( outrigger.Htot+1,  0, - outrigger.Hstory   )
    fig, ax = plt.subplots(figsize=(4, 5))   
    ax.plot(outrigger.Dw * 1000 , yVal ,  label='Coupled') 
    ax.plot(outrigger.DwC * 1000 , yVal,  label='Cantilever') 

    ax.set_title('Deflection')
    
    
    ax.set_xlabel ("Deflection")
    ax.set_ylabel ("Height")

    ax.grid (visible=True ,axis='both')
    ax.set_aspect ( math.ceil(max ( 1000*outrigger.DwC) / 10) * 10 / outrigger.Htot *2.5 )
    ax.set_ylim (0,  outrigger.Htot )
    ax.set_xlim (0,  math.ceil(max ( 1000*outrigger.DwC) / 10) * 10)
    ax.legend()
    # Save the plot to a BytesIO buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)  
    
    # Serve the image as a response
    return base64.b64encode(buf.read()).decode("utf-8")    


def Draw_Moment (outrigger):
    yVal = np.arange( outrigger.Htot+1,  0, - outrigger.Hstory   )
    xVal  = -outrigger.Mw
    fig, ax = plt.subplots(figsize=(4, 5))   
    ax.plot(xVal  , yVal ,  label='Coupled') 
    x1 = math.floor(min ( xVal) / 10000)*10000
    xVal  = -outrigger.MwC
    ax.plot(xVal , yVal, label='Cantilever') 

    # ax.set_aspect(1000)
    ax.set_aspect ( math.ceil(max ( xVal) / 10000) * 10000 / outrigger.Htot *2.5  )

    # ax.set_xlim ( math.floor(min ( xVal) /10000) * 10000, math.ceil(max ( xVal) /10000) * 10000)
    ax.set_xlim ( x1 , math.ceil(max ( xVal) /10000) * 10000)

    ax.set_ylim (0,  outrigger.Htot )
    ax.legend ()

    ax.set_title('Wall Moment')         
    
    ax.set_xlabel ("Moment")
    ax.set_ylabel ("Height")
    # ax.set_xticks  = [ 0   , 20000    , 30000 , 50000      ]
    ax.grid (visible=True ,axis='both')


    # Save the plot to a BytesIO buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)
    
    # Serve the image as a response
    return base64.b64encode(buf.read()).decode("utf-8")    


def Draw_Axial (outrigger):    
    yVal = np.arange( outrigger.Htot+1,  0, - outrigger.Hstory   )
    fig, ax = plt.subplots(figsize=(4, 5))  
    ax.plot(outrigger.Tc  , yVal) 
    # ax.plot(outrigger.TcC  , yVal) 
    ax.set_title('Axial Force')

    
    ax.set_xlim (0, math.ceil(max ( outrigger.Tc) /100) * 100)
    ax.set_ylim (0,  outrigger.Htot )
    ax.set_xlabel ("Axial Force")
    ax.set_ylabel ("Height")
    ax.grid (visible=True ,axis='both')
    ax.set_aspect ( math.ceil(max ( outrigger.Tc) /100) * 100 / outrigger.Htot *2.5 )
    ax.set_box_aspect( aspect=3)

    # Save the plot to a BytesIO buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)
    
    # Serve the image as a response
    return base64.b64encode(buf.read()).decode("utf-8")    
    # return Response(buf, mimetype='image/png')
    
if __name__ == '__main__':
     app.run(debug=True)