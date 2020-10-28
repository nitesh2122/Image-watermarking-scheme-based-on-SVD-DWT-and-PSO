% generating alpha using PSO

% velocity -  0 - 40
% ini vel - 20
% swarm range 1 - 10 
% inertia 0.5

alpha = PSO();
%alpha = 0.2;
%alpha = 0;

disp(alpha);
% image watermarkin
[SLL,ULLMA, VLLMA] = WaterMark(alpha);

% image extraction
Extracting_WaterMark(alpha,SLL,ULLMA,VLLMA);



