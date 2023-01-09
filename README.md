# PX913 Group Project
Group project with big man Jake

So basically, use this readme to update/explain work you've done so far. Feel free to play around with the structure if u like.

----------------------------------------------------------------------------------------------
------------------------------------Joseph:---------------------------------------------------
1. I made github rep :)




----------------------------------------------------------------------------------------------
-------------------------------------Jake:----------------------------------------------------
1. CURRENT POSITION: Have included;  DensityGrid, PhiGrid(DesnityGrid) w/ convergence criteria, EFieldGrid(PhiGrid).
2. Believe current issues relate to how we are initialising inner-values in PhiGrid.
3. Have emailed Heather to ask for info on initialization.

----------------------------------------------------------------------------------------------
1. I believe I've solved the issue. Not completely sure it's correct but we can compare to others.
2. Solved by setting N to be REAL64 (was previously INT) in phi grid. 
3. Also I enabled system to converge without use of ABS() in error quantifier (drms) by doing one run above the loop, and setting all values within to negative (line 166).
4. Haven't really bothered looking into why this works, so cant say for sure is a perfect system. 
