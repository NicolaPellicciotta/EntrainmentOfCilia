function bb= calibration_remove_nan(aa)
aa(aa ==0 | aa ==-1 )=nan;
bb=aa;
end