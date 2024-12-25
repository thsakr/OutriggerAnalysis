from flask import Flask, request, render_template , jsonify

import math     as mt
import numpy	as np

import datetime as dt



class   Outrigger:
    
    def __init__(self, data):
        self.sc =float ( data['sc'])
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

        MwC = np.zeros (Nst+1)
        TcC = np.zeros (Nst+1)
        DwC = np.zeros (Nst+1)
       
        for  iii in range (0,Nst+1,1 ):	
            xh = iii*self.Hstory
            z = xh/self.Htot
            MwC[iii] = -self.Load*(xh*xh)/2             
            DwC[iii] = self.Load * self.Htot**4 / (24 * EIcore )  * (3-4*z+ z**3)
          #  Dw[iii] = DwC[iii]

        Mw = np.zeros (Nst+1)
        Tc = np.zeros (Nst+1)
        Dw = np.zeros (Nst+1)

        FlO =np.zeros (n+1,dtype=int)
        self.Res    =   np.zeros((Nst+1,7),dtype=float)

        for  iii in range (0,n ):	
            FlO [iii] = int (self.x[iii]/self.Hstory)

        FlO[n]  =  Nst

        
        Macc = 0
        for iiii in range (0 , FlO[1]+1):                
                Mw[iiii] = MwC[iiii] 

        for  iii in range (0,n,1 ): 
            Macc = Macc + M[iii]                
            for iiii in range (FlO [iii] , FlO[iii+1]+1):                
                Mw[iiii] = MwC[iiii] + Macc 
                Tc[iiii] = TcC[iiii] + Macc / (self.D)


        for  iii in range (0,n,1 ):                         
            for iiii in range (0 , Nst+1 ):
                xh  =   iiii*self.Hstory
                xb  =   self.Htot-  xh
                z = xh/self.Htot
                if  (xh > self.x[iii] ):            # down
                    Dw[iiii]  = Dw[iiii] + M[iii] *xb**2 /(2*EIcore )
                   # print (xh,xb,Dw[iiii])
                else:                               # Up
                    Dw[iiii]  = Dw[iiii] +  M[iii] *(self.Htot -  self.x[iii])**2 /(2*EIcore )+ M[iii] *(self.Htot -  self.x[iii]) *(xb -(self.Htot-self.x[iii] ) )/EIcore 
                   # print (xh,xb,Dw[iiii])

        for iiii in range (0 , Nst+1 ):
            Dw[iiii] = DwC[iiii] - Dw[iiii]


            # output method 1
        self.Res = "<h3> Outrigger System Results </h3>"

        self.Res  = self.Res +   " <table  border='2' bgcolor='#efefef'> " 
        self.Res  = self.Res +   " <thead> "

        self.Res  = self.Res +   " <tr> " 
        
        self.Res  = self.Res +   " <th>  Height  </th> "

        self.Res  = self.Res +   " <th width = '150' >  Moment  </th> "
        self.Res  = self.Res +   " <th  width = '150'> Tension  </th> "
        self.Res  = self.Res +   " <th  width = '150'>  Deflection (mm) </th> "
        self.Res  = self.Res +   " </tr> "
        self.Res  = self.Res +   " </thead> "

        for ix in range (0,Nst+1):
            self.Res  = self.Res + " <tr> "
            self.Res  = self.Res + "<th> "+ str ("{:.2f}".format( self.Htot-  ix * self.Hstory )) +"  </th> "
            self.Res  = self.Res + "<th> "+ str (  "{:.2f}".format( Mw[ix] )) +"  </th> "
            self.Res  = self.Res + "<th> "+ str ( "{:.2f}".format( Tc[ix]) )+"  </th> "
            self.Res  = self.Res + "<th> "+str ( "{:.2f}".format( Dw[ix] * 1000 )) +"  </th> "                        
            self.Res  = self.Res + " </tr> "
            
        self.Res  =   self.Res   + "</table> "

        self.Res  =   self.Res   + "<h5> Cantilever Moment = " +str ( "{:.2f}".format( MwC[Nst] * 1000 ))+"</h5>"
        self.Res  =   self.Res   + "<h5> Cantilever Deflection = " +str ( "{:.2f}".format( DwC[0] * 1000 ))+" mm </h5>"
          

        
         # output method 2
        # self.res = {
        #         "Cantilever Moment": MwC,
        #         "Coupled Moment": Mw,
        #         "Cantilever Tension": TcC,
        #         "Coupled Tension": Tc,
        #         "Cantilever Deflection": DwC,
        #         "Coupled Deflection": Dw,
        #     }
        
            # output method 3
        # self.Res = "<h3> Outrigger System Results </h3>"

        # self.Res  = self.Res +   " <h4> Cantilever Moment </h4> " 
        # for ix in range (0,Nst+1):            
        #     self.Res  =   self.Res   + str (ix * self.Hstory ) +"  , "+ str (  MwC[ix] ) + "  ,  "
            
        # self.Res  =  self.Res   + " <h4> Coupled Moment </h4>  "     
        # for ix in range (0,Nst+1):            
        #     self.Res  =   self.Res   + str (ix * self.Hstory ) +"  , "+ str (  Mw[ix] ) + "  ,  "
        

        
        # self.Res  = self.Res+ " <h4> Cantilever Tension </h4>   " 
        # for ix in range (0,Nst+1):            
        #     self.Res  =   self.Res   + str (ix * self.Hstory ) +"  , "+ str (  TcC[ix] ) + "  ,  "
            
        # self.Res  =  self.Res   + "Coupled Tension  "     
        # for ix in range (0,Nst+1):            
        #     self.Res  =   self.Res   + str (ix * self.Hstory ) +"  , "+ str (  Tc[ix] ) + "  ,  "
        
          
        # self.Res  = self.Res+ " <h4> Cantilever Deflection </h4>  " 
        # for ix in range (0,Nst+1):            
        #     self.Res  =   self.Res   + str (ix * self.Hstory ) +"  , "+ str (  DwC[ix] ) + "  ,  "
            
        # self.Res  =  self.Res   + " <h4>  Coupled Deflection </h4>  "     
        # for ix in range (0,Nst+1):            
        #     self.Res  =   self.Res   + str (ix * self.Hstory ) +"  , "+ str (  Dw[ix] ) + "  ,  "
          


        return self.Res


app = Flask(__name__)

@app.route  ('/')

def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])


def analyze():
    My_Data = request.form.to_dict()
    #data = request.json
    outrigger = Outrigger(My_Data)
    result = outrigger.Out_Analysis()
    # return result

    # return jsonify({"status": "error", "message": "Analysis failed.", "details": result})

    return  result


if __name__ == '__main__':
     app.run(debug=True)