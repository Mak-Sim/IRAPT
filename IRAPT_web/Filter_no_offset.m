function [Frame]=Filter_no_offset(Imp,Frame,H_size,Ground)
Frame=Frame(:)';
Frame_ex=[Frame zeros(1,H_size)+Ground];
Frame=filter(Imp,1,Frame_ex);
Frame=Frame(H_size+1:end);
