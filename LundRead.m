function sig= LundRead(name, N)

if nargin< 2,
   error('sig= LundRead(name, N)');
end

% if name(1) ~= '_',
%    name= ['_',name];
% end
if name(end) ~= 'g',
   name= [name,'.ecg'];
end

fic= fopen(name, 'r');
sig= fread(fic, Inf, 'int16');
fclose(fic);

sig= sig(256+1:end);
sig= reshape(sig, [N, length(sig)/N])'; 
 
    