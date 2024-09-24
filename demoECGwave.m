
disp('Cluster Wave Analysis'); 
disp('Demo for Lund ECG and MIT ECG');

disp('Menu:');
disp('   1.- Lund case1 with 3 leads');
disp('   2.- Lund case1 with 9 leads');
disp('   3.- Lund case3 with 3 leads');
disp('   4.- Lund case33 with 3 leads');
disp('   5.- MIT database (100.hea) with 2 lead');
disp('   6.- MIT database (208.hea) with 2 lead');
cas= input('Input:');
switch cas
    case 1
       sig= LundRead('./Lundcases/case1_3L.ecg',3);
    case 2
       sig= LundRead('./Lundcases/case1_9L.ecg',9); 
    case 3
       sig= LundRead('./Lundcases/case3_3L.ecg',3); 
    case 4
       sig= LundRead('./Lundcases/case33_3L.ecg',3); 
    case 5
       hea= sopen('./MITcases/100.hea');
       sig= sread(hea);
    case 6
       hea= sopen('./MITcases/208.hea');
       sig= sread(hea);
end
L= floor(size(sig,1)/4);
sig= sig(1:L, :);
ECGwave(sig);