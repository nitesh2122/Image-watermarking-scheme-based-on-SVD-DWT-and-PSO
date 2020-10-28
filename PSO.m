function alpha = PSO()
    disp('func - PSO')
    swarm = [0.1,3,1,7,9];
    vel = [20,20,20,20,20];
    Iw = 0.5;

    pbest = swarm;
    pfitness = [-1000000,-1000000,-1000000,-1000000,-1000000];
    gbest = 1;
    
    ind = 0;
    for j=1:2
        for i=1:5 
            ind = ind + 1;
            disp(ind);
            temp = fit(swarm(i));
            if( temp > pfitness(i) )
                pbest(i) = swarm(i);
                pfitness(i) = temp;
            end
        end
       
        for i=1:5
            if( pfitness(i) > pfitness(gbest) )
                gbest = i;
            end
        end
        
        for i=1:5
            % finding velocity
            % updating position
            vel(i) =  Iw * vel(i) + 2 * rand(1,1) * (pbest(i)-swarm(i)) + 2*rand(1,1)*(pbest(gbest) - swarm(i)) ;
            if vel(i)>40
                vel(i)=0;
            end
            if vel(i)<0
                vel(i)=0;
            end
            swarm(i) = swarm(i) + vel(i);
            
        end
    end
    
    alpha = pbest(gbest);
end


function fitness=fit(alpha)
    disp('func - PSO_fit')
    image = imread('Image.bmp'); image = imresize(image,[1024,1024]);
    watermark = imread('Watermark.bmp'); watermark = imresize(watermark,[128,128]);
    
    [ULLMA,VLLMA,SLL,watermarked] = watermarker(alpha);
    extracted = extractor(alpha, ULLMA, VLLMA, SLL,watermarked);
    
    a = normxcorr2(image,watermarked);
    b = normxcorr2(watermark,extracted);
    
    fitness = max(a(:)) + max(b(:));
end


function [x,y,z,newImage] = watermarker(alpha)
    disp('func - PSO_Watermark')
    image = imread('Image.bmp'); Image = imresize(image,[1024,1024]);
    watermark = imread('Watermark.bmp'); WaterMark = imresize(watermark,[128,128]);
    
    % dwt over image
    [LL, LH, HL, HH] = dwt2(Image, 'haar') ;
    
    % divide LL into 4 * 4 non overlapping blocks denote by and apply svd on them 
    % apply svd
    
    ULL=zeros(4,4,128*128);
    SLL=zeros(4,4,128*128);
    VLL=zeros(4,4,128*128);
    k = 1;
   
    for i = 1:4:512
        for j = 1:4:512
            mat=zeros(4,4);
            for ii = 1:4
                for jj = 1:4
                    mat(ii,jj) = LL(i+ii-1,j+jj-1);
                end
            end
            [ULL(:,:,k),SLL(:,:,k),VLL(:,:,k)] = svd(mat); k=k+1;
        end
    end
        
    % modify the first and largest singular value of each block by a pixel
    % of the watermarked image 
    % SLLM[i,j] = SLL[i,j](1,1) + alpha * w[i,j];
    SLLM=zeros(4,4,128*128);
    for i = 1:128*128
        dummy = SLL(:,:,i);
        ii = floor((i-1)/128)+1;
        jj = mod((i-1),128)+1;
        abc = alpha * WaterMark(ii, jj);
        dummy(1,1) = double(dummy(1,1)) + double(abc);
        SLLM(:,:,i) = dummy;
    end
    
    % apply block svd to each block SLLM[i,j] 
    ULLMA = zeros(4,4,128*128);
    SLLMA = zeros(4,4,128*128);
    VLLMA = zeros(4,4,128*128);
    for i = 1:128*128 
        [ULLMA(:,:,i),SLLMA(:,:,i),VLLMA(:,:,i)] = svd(SLLM(:,:,i));
    end
    
    % apply inverse block svd on SLLMA using orthogonal matric ULL, VLL
    subLLM = zeros(4,4,128*128);
    for i = 1:128*128 
        temp = VLL(:,:,i)';
        subLLM(:,:,i)= ULL(:,:,i) * SLLMA(:,:,i) * temp ;
    end

    LLM = zeros(512,512);

    jj = 1;
    for i = 1:128*128 
        ii = 1 + (floor((i-1)/128))*4;
        for iii= 0:3
            for jjj = 0:3
                   LLM(ii+ iii,jj+jjj) = subLLM(iii+1,jjj+1,i);
            end
        end
        jj = jj + 4 ;
        jj = mod(jj,512);
    end
    
      
    % apply inverse dwt on the modified LLM with other subbands
    newImage = uint8(idwt2(LLM, LH, HL, HH, 'haar'));
    x = ULLMA;
    y = VLLMA;
    z = SLL;
end

function extracted = extractor(alpha,ULLMA,VLLMA,SLL, watermarked)

    disp('func - PSO_extractor');
    
    [LL,LH,HL,HH] = dwt2(watermarked, 'haar');
    
    % divide LL into 4 * 4 non overlapping blocks denote by and apply svd on them 
    % apply svd
    
    ULL = zeros(4,4,128*128);
    SLLMA = zeros(4,4,128*128);
    VLL = zeros(4,4,128*128);
    
    k = 1;
    for i = 1:4:512 
        for j = 1:4:512
            mat = zeros(4,4);
            for ii = 1:4
                for jj = 1:4
                    mat(ii,jj) = LL(i+ii-1,j+jj-1);
                end
            end
            [ULL(:,:,k),SLLMA(:,:,k), VLL(:,:,k)]  = svd(mat);  k=k+1;
        end
    end
    
    
    SLLM = zeros(4,4,128*128);
    for i = 1:128*128
        SLLM(:,:,i) = ULLMA(:,:,i) * SLLMA(:,:,i)* (VLLMA(:,:,i)');
    end
    
    w = zeros(128,128);
    for i=1:128*128 
        a = floor((i-1)/128)+1;
        b = mod((i-1),128)+1;
        p = double(SLLM(1,1,i));
        q = double(SLL(1,1,i));
        w(a,b) = (p - q) / alpha ;
    end
    
    extracted = uint8(w);
end