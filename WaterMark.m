function [w,x,y] = WaterMark(alpha) 
    disp('func - watermark');
    % read image 
    image = imread('Image.bmp');
    image = imresize(image, [1024,1024]);

    % read watermark
    WaterMark = imread('Watermark.bmp');
    WaterMark = imresize(WaterMark,[128,128]);
    
    
    % dwt over image
    [LL, LH, HL, HH] = dwt2(image, 'haar');
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
    %figure(1);imshow(uint8(LLM));
    
    % apply inverse dwt on the modified LLM with other subbands
    WaterMarked = idwt2(LLM, LH, HL, HH, 'haar');
    %imshow(uint8(WaterMarked));
    imwrite(uint8(WaterMarked),'watermarkedImage.bmp');
    
    disp(psnr(image,uint8(WaterMarked)));
    robust = normxcorr2(image,WaterMarked);
    disp(max(robust(:)));
    w = SLL;
    x = ULLMA;
    y = VLLMA;

end