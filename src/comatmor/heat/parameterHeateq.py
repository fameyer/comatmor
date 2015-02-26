# mat-File name for matrix input
# make $HOME available !!!
# Arbitrary parameterranges possible - Keep attention to correct definitions!
matfile = {'Kc1': ('/home/310191226/pymorDir/comatmor/src/comatmor/heat/Kc1.mat',
['c1'],[1],[(1,50)]),'Kc2': ('/home/310191226/pymorDir/comatmor/src/comatmor/heat/Kc2.mat',['c2'],[1],[(1,50)]),'Lc1': ('/home/310191226/pymorDir/comatmor/src/comatmor/heat/Lc1.mat',['c1'],[1],[(1,50)]),'Lc2': ('/home/310191226/pymorDir/comatmor/src/comatmor/heat/Lc2.mat',['c2'],[1],[(1,50)])}
stiffNames = ('Kc1','Kc2')
rhsNames = ('Lc1','Lc2') 
u0file = {'u0': '/home/310191226/pymorDir/comatmor/src/comatmor/heat/u0.mat'}
massfile = {'Dc': '/home/310191226/pymorDir/comatmor/src/comatmor/heat/Dc.mat'}
trainingSetfile = {'training_set': '/home/310191226/pymorDir/comatmor/src/comatmor/heat/training_set.mat'}
