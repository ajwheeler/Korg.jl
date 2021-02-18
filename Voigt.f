      real*8 function voigt(a,v)
c******************************************************************************
c     This routine calculates the voigt function using the approximations   
c     found in the *new series-astronomy and astrophysics* volume of        
c     landolt-bornstein.    
c******************************************************************************
      implicit real*8 (a-h,o-z) 

      a2 = a*a
      v2 = v*v
      if(a .eq. 0.) go to 62 
      if(a .le. 0.2) go to 50
      if(a.le.1.4 .and. a+v.le.3.2) go to 30    


c  case 1, a .gt. 1.4 or (a .gt. 0.2 and a+v .gt. 3.2)                   
      u = 1.4142136*(a2 + v2)                  
      voigt = 0.7978847*a/u*(1.0 + (3.0*v2-a2)/u**2 + (15.0*v2*v2      
     1 -30.0*a2*v2+3.0*a2*a2)/u**4)                                   
      voigt = voigt/1.772454                                         
      return                                                        
30    u=0.979895023-0.962846325*a+0.532770573*a2-0.122727278*a*a2  
31    h0 = dexp(-v2)                                              
      if(v .ge. 2.4) go to 35                                    
      if(v .ge. 1.3) go to 33                                   


c  case 2,  0.2 .lt. a .le. 1.4 and a+v .le. 3.2                         
      h1=-1.12470432-0.15516677*v+3.28867591*v2-2.34357915*v*v2
     1 +0.42139162*v2*v2                                      
      go to 36                                               
33    h1=-4.48480194+9.39456063*v-6.61487486*v2+1.98919585*v*v2        
     1 -0.2204165*v2*v2                                               
      go to 36                                                       
35    h1=(0.554153432+0.278711796*v-0.188325687*v2+0.042991293*v*v2
     1 -0.003278278*v2*v2)/(v2 - 1.5)                             
36    h2 = (1.0 - 2.0*v2)*h0                                     
      if(a .le. 0.2) go to 52                                   
      h1 = h1 + 1.1283790*h0                                   
      h2p = h2                                                
      h2 = h2 - h0 + 1.1283790*h1                            
      h3 = 0.37612635*(1.0-h2p) - 0.6666667*v2*h1 + 1.1283790*h2       
      h4 = 0.6666667*v2*v2*h0 - 0.37612635*h1 + 1.1283790*h3          
      voigt = u*(h0 + h1*a + h2*a2 + h3*a*a2 + h4*a2*a2)             
      voigt = voigt/1.772454                                        
      return                                                       


c  case 3,  a .le. 0.2 and v .lt. 5.0                                    
50    if(v .ge. 5.0) go to 60                                     
      go to 31                                                   
52    voigt = h0 + h1*a + h2*a2                                 
      voigt = voigt/1.772454                                   
      return                                                  


c  case 4,  a .le. 0.2 and v .ge. 5.0                                    
60    voigt = a/(1.772454*v2)*(1.0 + 1.5/v2 + 3.75/(v2*v2))  
      voigt = voigt/1.772454                                
      return                                               


c  case 5, a .eq. 0.0                                                    
62    voigt = exp(-v2)/1.772454                           


      return                                             
      end                                               
