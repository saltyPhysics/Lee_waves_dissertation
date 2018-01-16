# Lee waves dissertation

Here are a collection of python notebooks and some modules they utilize for internal lee wave analysis of a 100km transect
in the Southern Ocean past an island pack Shag Rocks. 

## **Contents**
* Internal wave analysis at Shag Rocks - https://github.com/saltyPhysics/Lee_waves_dissertation/blob/master/Internal_waves_shag_rocks.ipynb

* Potential Energy for Internal Waves - https://github.com/saltyPhysics/Lee_waves_dissertation/blob/master/Potential-Energy-For-Internal-Waves.ipynb

* Internal Wave Calculation functions - https://github.com/saltyPhysics/Lee_waves_dissertation/blob/master/internal_waves_calculations.py

* General Oceanographic Tools - https://github.com/saltyPhysics/Lee_waves_dissertation/blob/master/oceans.py


![img](http://latex.codecogs.com/svg.latex?%2A%2ASETUP%2A%2A%0D%0A%0D%0A%23%23%23%23+x%2C+y%2C+and+z+steps%3A%0D%0A%0D%0A%24%24+%5CDelta+x++%3D+%5CDelta+t+%5Cbigg+%5B+%5Cbigg+%28%5Cfrac%7BN%5E2+k%7D%7Bm%5E2+%5COmega%7D+%5Cbigg+%29++%2B+U%28x%2Cy%2Cz%2Ct%29+%5Cbigg+%5D+%24%24%0D%0A%0D%0A%24%24+%5CDelta+y++%3D+%5CDelta+t+%5Cbigg+%5B+%5Cbigg+%28+%5Cfrac%7BN%5E2+l%7D%7Bm%5E2+%5COmega%7D+%5Cbigg+%29+%2B+V%28x%2Cy%2Cz%2Ct%29+%5Cbigg+%5D+%24%24%0D%0A%0D%0A%24%24+%5CDelta+z++%3D+%5CDelta+t+%5Cbigg+%5B+%5Cfrac%7BN%5E2%28k%5E2+%2B+l%5E2%29%7D%7Bm%5E3+%5COmega%7D+++%5Cbigg+%5D+%24%24%0D%0A%0D%0A%23%23%23%23+Get+change+in+position+and+new+field+values%3A+_Inverse+Haversine+Forumala+for+spherical+distance_%0D%0A%0D%0A%24%24+%5CDelta+longitude++%3D+%28%5CDelta+x%29%28%2A%2Aconversion+factor%2A%2A%29++%24%24%0D%0A%0D%0A%24%24+%5CDelta+latitude++%3D+%28%5CDelta+y%29%28111.11+km%29++%24%24%0D%0A%0D%0A%23%23%23%23+Changes+in+wave+parameters%3A%0D%0A%0D%0A%0D%0A%0D%0A%24%24+%5CDelta+k++%3D+++%5CDelta+t+%5Cbigg+%5B++r_x++%2B+k%5Cfrac%7B%5Cpartial+U%7D%7B%5Cpartial+x%7D+%2B+l%5Cfrac%7B%5Cpartial+V%7D%7B%5Cpartial+x%7D+%5Cbigg+%5D+%24%24%0D%0A%0D%0A%24%24+%5CDelta+l++%3D+++%5CDelta+t+%5Cbigg+%5B++r_y++%2B+k%5Cfrac%7B%5Cpartial+U%7D%7B%5Cpartial+y%7D+%2B+l%5Cfrac%7B%5Cpartial+V%7D%7B%5Cpartial+y%7D+%5Cbigg+%5D+%24%24%0D%0A%0D%0A%24%24+%5CDelta+m++%3D+++%5CDelta+t+%5Cbigg+%5B++r_z++%2B+k%5Cfrac%7B%5Cpartial+U%7D%7B%5Cpartial+z%7D+%2B+l%5Cfrac%7B%5Cpartial+V%7D%7B%5Cpartial+z%7D+%5Cbigg+%5D+%24%24%0D%0A%0D%0A%0D%0Awhere+%24r_%7Bx%2Cy%2Cz%7D%24+is+the+refraction+index%3A%0D%0A%0D%0A%24%24+r+%3D+%5Cfrac%7BN%28k%5E2+%2Bl%5E2%29%7D%7Bm%5E2+%5COmega%7D+%5Cfrac%7B%5Cpartial+N%7D%7B%5Cpartial+%28x%2Cy%2Cz%29%7D+%24%24%0D%0A%0D%0A%23%23%23%23+Wave+dispersion+relation+change%0D%0A%0D%0AChanges+to+intrinsic+frequency%2Fdispersion+relation+along+the+wave+ray+is+equal+to+the+sum+of+the+changes+to+the+wave+number+components+%3A%0D%0A%0D%0A%24%24+%5CDelta+%5COmega+%3D+%5CDelta+t+%5Cbigg+%5B+%28r_x+%2B+r_y+%2B+r_z%29+%2B+k%5CDelta+U+%2B+l%5CDelta+V+%5Cbigg+%5D+++%24%24)
