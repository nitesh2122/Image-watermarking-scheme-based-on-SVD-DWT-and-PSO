function Extracting_WaterMark(alpha, SLL, ULLMA, VLLMA)
    
    disp('func - Extracting_watermark')
    WaterMarked = imread('watermarkedImage.bmp');
    [LL,LH,HL,HH] = dwt2(WaterMarked, 'haar');
    
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
    
    % write w
    imwrite(uint8(w),'extractedWatermark.bmp');
end