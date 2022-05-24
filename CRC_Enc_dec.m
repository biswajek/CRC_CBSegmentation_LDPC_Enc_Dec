%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %BISWAJEET KUMAR%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Generate TB of size 20496
msg = randi([0,1],1,20496);

%CRC Polynomial 1 %
crc_poly1=[1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
%CRC Polynomial 2 %
crc_poly2=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];

% Transport Block CRC calculation %
[M,N] = size(crc_poly1);
msg = [msg,zeros(1,N-1)];

%Append the code block CRC
[q,r] = deconv(msg,crc_poly1);

r=abs(r);

for i=1:length(r)
    a=r(i);
    if ( mod(a,2)== 0 )
        r(i)=0;
    else
        r(i)=1;
    end
end

%Attached TB CRC to get the code block
codeBlock = bitor(msg,r);

[A,B] = size(codeBlock);

%Base graph selection %
base_graph = 'BG1';
%CRCB Length
L = 24;

% Code Block Segmentation %
if base_graph == 'BG1'
    Kcb = 8448; % Max code block size for base graph 1
    C = ceil(B/(Kcb - L));
    Bdash = B + C * L;
    Kdash = Bdash/C;
    % Kb for Base Graph 1
    Kb = 22;
else
    Kcb = 3840; % Max code block size for base graph 2
    C = ceil(B/(Kcb - L));
    Bdash_BG2 = B + C * L;
    Kdash = Bdash/C;
    %Kb for Base Graph 2
    if B > 640
        Kb = 10;
    elseif B > 560
        Kb = 9;
    elseif B > 192
        Kb = 8;
    else
        Kb = 6;
    end
end

% Z0 = [2 4 8 16 32 64 128 256];
% Z1 = [3 6 12 24 48 96 192 384];
% Z2 = [5 10 20 40 80 160 320];
% Z3 = [7 14 28 56 112 224];
% Z4 = [9 18 36 72 144 288];
% Z5 = [11 22 44 88 176 352];
% Z6 = [13 26 52 104 208];
% Z7 = [15 30 60 120 240];
 
Z = [2 4 8 16 32 64 128 256,3 6 12 24 48 96 192 384,5 10 20 40 80 160 320,7 14 28 56 112 224,...
    9 18 36 72 144 288,11 22 44 88 176 352,13 26 52 104 208,15 30 60 120 240];
M = Kb.*Z;
I = find(M > Kdash); 
Zc = min(Z(I));

if base_graph == 'BG1'
    K = 22 * Zc;
else
    K = 10 * Zc;
end

cb = zeros(C,Kdash-L-1);
s = 1;
for c = 1:C
    for k = 1:Kdash-L
        cb(c,k) = codeBlock(s);
        s = s + 1;
    end
    
    if C > 1
        %Find the code block CRC using CRC24(B) polynomial
        [M,N] = size(crc_poly2);
        cb_temp = [cb(c,1:end),zeros(1,N-1)];
        %Append the code block CRC
        [q,r] = deconv(cb_temp,crc_poly1);
        r=abs(r);
        %finding the absolute value of remainder
        for i=1:length(r)
            a=r(i);
            if ( mod(a,2)== 0 )
                r(i)=0;
            else
                r(i)=1;
            end
        end
        cb_temp = bitor(cb_temp,r);
        ind = 1;
        for k=Kdash - L+1:Kdash
            cb(c,k) = cb_temp(ind);
            ind = ind + 1;
        end
    end
    for k=Kdash+1: K
        cb(c,k) = 0;
    end
    
    nbg = 1; % Base graph 1
    nldpcdecits = 25; % Decode with maximum no of iteration

    %----------------------------------LDPC encoding----------------------- ------------------------

    ldpc_coded_bits = double(LDPCEncode(cb(c,1:end)',nbg)); 
    mod_output = 2*(ldpc_coded_bits-0.5);
    noise_power = (10^-5);
    noise = sqrt(noise_power)*randn(size(mod_output));
    rx_sig = mod_output + noise; % to add noise
    llr0 =  abs(-1 + rx_sig);   % BPSK demod
    llr1 =  abs(1 + rx_sig);    % BPSK demod
    llr = log(llr0./llr1);      % ldpc decoder requires log(p(r/0)/p(r/1))
    demod_output = llr;

    %----------------------------------LDPC decoding----------------------- -------------------------
    outputbits = double(LDPCDecode(demod_output,nbg,nldpcdecits));
    temp = cb(c,1:end)';
    errors = find(outputbits - temp)
end
