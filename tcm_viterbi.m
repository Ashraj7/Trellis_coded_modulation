clear all;
close all;
clc;
bits=[1 0 0 1 0 1 0 1 0 0 1 0
          1 1 0 0 1 1 0 1 0 1 1 0
          0 1 1 0 0 1 1 1 0 0 1 0
          1 0 1 1 0 0 1 0 1 1 1 0
          0 0 1 1 0 1 0 0 1 1 0 1
          1 1 0 0 1 0 1 1 0 0 1 0
          1 0 1 0 1 1 0 1 0 0 1 1
          0 1 0 1 1 0 0 1 1 0 1 0
          1 1 0 0 1 1 0 1 0 0 1 1
          0 1 0 1 1 0 1 0 0 1 1 0
          1 0 1 1 0 0 1 0 1 1 0 1
          0 1 1 0 0 1 1 0 1 0 1 1 ];
bits = reshape(bits.',1,[]);
numtribits = length(bits)/3;
tribits = reshape(bits,3,[])';
disp("Input bits:");
disp(bits);
disp("Grouped Tribits:");
disp(tribits);

tribit_values = zeros(1,numtribits);
for r = 1:numtribits
    value = tribits(r,1)*4 + ...
            tribits(r,2)*2 + ...
            tribits(r,3);
    tribit_values(r) = value;
endfor
tribit_values = [tribit_values 0];    %adding flushing tribit
disp("Tribit Values:");
disp(tribit_values);
fsm = [
0 8 4 12 2 10 6 14;
4 12 2 10 6 14 0 8;
1 9 5 13 3 11 7 15;
5 13 3 11 7 15 1 9;
3 11 7 15 1 9 5 13;
7 15 1 9 5 13 3 11;
2 10 6 14 0 8 4 12;
6 14 0 8 4 12 2 10];
state = 0;   %initial state is zero
constellations = zeros(1,length(tribit_values));
for i = 1:length(tribit_values)
    constellation = fsm(state+1, tribit_values(i)+1);
    constellations(i) = constellation;
    state = tribit_values(i);
endfor
disp("Constellation Points:");
disp(constellations);
constellation_map = [
 +1 -1;
 -1 -1;
 +3 -3;
 -3 -3;
 -3 -1;
 +3 -1;
 -1 -3;
 +1 -3;
 -3 +3;
 +3 +3;
 -1 +1;
 +1 +1;
 +1 +3;
 -1 +3;
 +3 +1;
 -3 +1];

mapped_symbols = [];

for i = 1:length(constellations)
    point = constellations(i);
    pair = constellation_map(point+1,:);
    mapped_symbols = [mapped_symbols pair];
endfor
disp("Mapped Symbols:");
disp(mapped_symbols);
interleave_table = [0;
                                   1;
                                   8;
                                   9;
                                   16;
                                   17;
                                   24;
                                   25;
                                   32;
                                   33;
                                   40;
                                   41;
                                   48;
                                   49;
                                   56;
                                   57;
                                   64;
                                   65;
                                   72;
                                   73;
                                   80;
                                   81;
                                   88;
                                   89;
                                   96;
                                   97;
                                   2;
                                   3;
                                   10;
                                   11;
                                   18;
                                   19;
                                   26;
                                   27;
                                   34;
                                   35;
                                  42;
                                  43;
                                  50;
                                  51;
                                  58;
                                  59;
                                  66;
                                  67;
                                  74;
                                  75;
                                  82;
                                  83;
                                  90;
                                  91;
                                  4;
                                  5;
                                  12;
                                  13;
                                  20;
                                  21;
                                  28;
                                  29;
                                  36;
                                  37;
                                  44;
                                  45;
                                  52;
                                  53;
                                  60;
                                  61;
                                  68;
                                  69;
                                  76;
                                  77;
                                  84;
                                  85;
                                  92;
                                  93;
                                  6;
                                  7;
                                  14;
                                  15;
                                  22;
                                  23;
                                  30;
                                  31;
                                  38;
                                  39;
                                  46;
                                  47;
                                  54;
                                  55;
                                  62;
                                  63;
                                  70;
                                  71;
                                  78;
                                  79;
                                  86;
                                  87;
                                  94;
                                  95];

interleaved_symbols = zeros(1,98);

for i = 1:98
    interleaved_symbols(i) =mapped_symbols(interleave_table(i)+1);
endfor
disp("Interleaved Symbols:");
disp(interleaved_symbols);

%%adding awgn
snr=10;
received_symbols = awgn(interleaved_symbols, snr, 'measured');
%%deinterleaving

deinterleaved_symbols = zeros(1,length(received_symbols));
for i=1:98
  deinterleaved_symbols(interleave_table(i)+1)= received_symbols(i);
endfor
disp("Deinterleaved Symbols:");
disp(deinterleaved_symbols);


%%viterbi
history=zeros(8,49);
survivor_tribit=zeros(8,49);
path_metric=[0 inf inf inf inf inf inf inf];
idx=1;
for j=1:49
  I_rx=deinterleaved_symbols(idx);
  Q_rx=deinterleaved_symbols(idx+1);
  temp_path=[inf inf inf inf inf inf inf inf];
  for state = 0:7
    if path_metric(state+1)==inf
      continue;
    endif
    for inp = 0:7
      next_state=inp;
      exp_const_pt=fsm(state+1,inp+1);
      exp_sym=constellation_map(exp_const_pt+1,:);
      I_exp= exp_sym(1);
      Q_exp=exp_sym(2);
      branch_metric =  (I_exp-I_rx)^2 + (Q_exp-Q_rx)^2;
      candidate_metric= path_metric(state+1)+branch_metric;
      if candidate_metric<temp_path(next_state+1)
        temp_path(next_state+1)=candidate_metric;
        history(next_state+1,j)=state;
        survivor_tribit(next_state+1,j)=inp;
      endif
    endfor
  endfor
  idx=idx+2;
  path_metric=temp_path;
endfor

[path,index]=min(path_metric);
current_state=index-1;
recovered_tribit=zeros(1,49);
for i= 49:-1:1
   recovered_tribit(i)= survivor_tribit(current_state+1,i);
   prev_state= history(current_state+1,i);
   current_state=prev_state;
 endfor
recovered_tribit = recovered_tribit(1:48);  %%remove the flushing tribit
disp("Recovered tribit values:");
disp(recovered_tribit);

for k = 1:48
    recov_bits(k,:) = dec2bin(recovered_tribit(k),3) - '0';
endfor
decoded_bits = reshape(recov_bits.',1,[]);
disp("Recovered bits:");
disp(decoded_bits);
bit_errors = sum(bits ~= decoded_bits);
disp("Bit errors:");
disp(bit_errors);
