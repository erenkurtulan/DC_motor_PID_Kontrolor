clc
clear
E=24;
Ra= 0.4;
La=2.7;
Ja=0.0004;
Ba=0.0022;
K=0.0015;           %%parametrelerin girişi
Kb=0.005;
R=1;  %%referans değer
Kt=10;
Kh=0.05;

%%isterler - geçici rejim kriterleri
ts_yerlesme=2; %%5to
istenen_asim=0.0432;
%%
s = tf('s');   %s değişkenini tanımla
n=K;
d= [(Ja*La) (Ja * Ra + La * Ba) (Ra * Ba + K *Kb)];
actf=tf(n,d);  %%açık çevrim tf bul
kok=roots(d);  %%kök bul    
tomin=1/max(abs(kok));  %% to seç
Ts=tomin/Kt;          %% örnekleme zamanı seç
%%
ayrikzaman=c2d(actf,Ts); %%sürekli zamandan ayrık zamana geçiş
[n_z,d_z]=tfdata(ayrikzaman,'v'); %% ayrık zaman tf fonk. vektörel yazdır
syms z;
Gz = poly2sym(n_z, z) / poly2sym(d_z, z);  % Transfer fonksiyonunu sembolik yap
ayrikzamankok=roots(d_z);
%%


%

%ileriyol=actf; %ileri yol fonk tanımla

tf_sistem= feedback(actf, 1); %sistemin kapalı çevrim tf bul

es_step= R/(1+actf); %basamak giriş için ess tanımla

ess_step=dcgain(es_step);

es_ramp= R/(s+ s * actf);

ess_ramp=dcgain(es_ramp);

fprintf('Basamak Giriş için Kararlı Durum Hatası: %f\n', ess_step);
fprintf('Rampa Giriş için Kararlı Durum Hatası: %f\n', ess_ramp);
%%


%%%% hesaplamalar
ksi= sqrt((log(istenen_asim)^2) / (pi^2 + log(istenen_asim)^2));%%ksi hesapla aşım formülünden
wn=4/(ksi*ts_yerlesme);
ess=2*ksi*Kh/wn;




%%%%%%%%%% istenilen 2. dereceden denklem

istenilen_denk= wn^2/(s^2+2*ksi*wn*s+wn^2);

istenen_gcz= c2d(istenilen_denk, Ts);




sigma=-4/ts_yerlesme; %benzetilecek kokun reel kısmını bul
teta=acos(ksi); %% reel kökle imajinerin açısını bul
imajiner_kok= sigma*1i * tan(teta); %% imajiner kökü bul

s1_kok=sigma+imajiner_kok; %%benzetilecek kökleri yaz
s2_kok=sigma-imajiner_kok;

z2_kok = exp(s1_kok * Ts); %% istenen s köklerinden istenen z köklerine geç
z1_kok= exp(s2_kok * Ts);
%% kd ve kp için gerekli parametreler
zkok_genlik=abs(z1_kok); %% istenen z kökünün genliğini bul
beta= atan2(imag(z1_kok), real(z1_kok));
Gz_kok1 = polyval(n_z, z1_kok) / polyval(d_z, z1_kok); % g(z1_kok) hesapla
Gz1_genlik = abs(Gz_kok1);  % Genliği al
Gz_kok1aci=angle(Gz_kok1);
kv=1/ess; %%hız hatası bul
%Ki=kv*Ts/(limit( Gz, z, 1)); %%ki hesapla
a=limit( Gz, z, 1);
Ki=kv*Ts/a;

Kd= (zkok_genlik/sin(beta))*(Ki*sin(beta)/(zkok_genlik-2*cos(beta)+(1/zkok_genlik))+sin(Gz_kok1aci)/Gz1_genlik);

Kp=-cos(Gz_kok1aci)/Gz1_genlik-(2*Ki*zkok_genlik)* (zkok_genlik-cos(beta))/(zkok_genlik^2-2*zkok_genlik*cos(beta)+1)+(sin(Gz_kok1aci)*-zkok_genlik+cos(beta)*sin(Gz_kok1aci))/(Gz1_genlik)*sin(beta);
%Kp= - (cos(Gz_kok1aci)/Gz1_genlik) - ( ( 2*Ki*zkok_genlik ) * ( zkok_genlik-cos(beta) ) / ( zkok_genlik^2- 2*zkok_genlik*cos(beta)+1 ) ) +  ( ( -zkok_genlik*sin(Gz_kok1aci)+cos(beta)*sin(Gz_kok1aci) ) / ( Gz1_genlik*sin(beta) ) );
Ki = double(Ki); 
Kd = double(Kd); 
Kp = double(Kp); 
fprintf('G(z1) acisi: %f\n', Gz_kok1aci);
fprintf('G(z1) genlik: %f\n', Gz1_genlik);
fprintf('beta: %f\n', beta);
fprintf('z1 kok genlik: %f\n', zkok_genlik);
fprintf('kd: %f\n', Kd);
fprintf('kp: %f\n', Kp);
fprintf('ki: %f\n', Ki);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

