function Iqrs= QRSdect(sig)
% function Iqrs= QRSdect(sig)

delta=5;
X= [sig(1:end-delta), sig(1:end-delta)-sig(delta+1:end)];
msig= mean(X(:,1));
X(:,1)= X(:,1)-msig;
Iqrs= find(X(:,1)>0);
if Iqrs< length(sig)/2, Iqrs= find(X(:,1)<0); end
figure; plot(X(:,1),X(:,2));
hold on; plot(X(Iqrs,1), X(Iqrs,2), '.r');
figure; plot(sig)
hold on; plot(sig(Iqrs),'.r');


