load X1.txt;
load Y1.txt;
load T1.txt;

clf;
nf=40;

fp=fopen('ontheway.bin');
sz=fread(fp,[1,3],'int32'); L=sz(3);M=sz(1);D=sz(2);
T =reshape(fread(fp,L*M*D,'double'),sz);
fclose(fp);

for l=1:nf
  plot(X1(:,1),  X1(:,2),  'bo','MarkerSize',8); hold on;
  plot(T (:,1,l),T (:,2,l),'ro','MarkerSize',5,'MarkerFaceColor',[1,0,0]);
  axis([-2.5,2.5,-2.5,2.5]);
  fn=sprintf('ontheway-%05d.png',l);
  print(fn,'-dpng');
  clf;
end

%%axis([-4.5,4.5,-1.5,1.5]);
