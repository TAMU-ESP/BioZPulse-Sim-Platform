function RGB2=rgb_edit_fn(RGB,A)
R = RGB(:,:,1); % turn your data into a "pseudo-gray" rgb image.
G = RGB(:,:,2);
B = RGB(:,:,3);
R(isnan(A)) = 1; %turn NaNs(I) to Red
G(isnan(A)) = 0;
B(isnan(A)) = 0;

R(A==inf) = 0; %turn infs(V) to Black
G(A==inf) = 0;
B(A==inf) = 0;

R(A==-inf) = 0; %turn infs(V) to Black
G(A==-inf) = 1;
B(A==-inf) = 0;

RGB2(:,:,1) = R; %combine color slices into A
RGB2(:,:,2) = G;
RGB2(:,:,3) = B;

RGB2=RGB;