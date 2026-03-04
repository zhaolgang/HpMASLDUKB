# HpMASLDUKB
Code for Associations between H pylori and MASLD risk in the UK Biobank

/*------------------------------------------------------------------------------------------------------------------------------*/
/*       Project name: HPYUKB                                                                                                   */
/*          Objective: To explore Helicobacter pylori and liver outcomes in the UKB                                             */
/*            Dataset: Hpylori, covariates, medication, outcomes, biomarkers                                                    */
/*         Programmer: Longgang Zhao (nhlon@channing.harvard.edu)                                                               */
/*   Program reviewer: Xinyuan Zhang (hpxzh@channing.harvard.edu)                                                               */
/*               Date: 12/08/2023 -- 01/30/2024                                                                                 */
/*           Exposure: Helicobacter pylori                                                                                      */
/*            Outcome: NAFLD, cirrhosis, severe liver disease, and liver death                                                  */
/*         Covariates: Age, sex, ethnicity, TDI, body mass index, diabetes, physical activity                                   */
/*------------------------------------------------------------------------------------------------------------------------------*/

/*Macros and notes*/
%macro histogram (datain,varlist);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne );%let var=%scan(&varlist.,&varnumber.);
title "Distribution of &var";proc sgplot data=&datain; histogram &var;run;
%let varnumber=%eval(&varnumber.+1); %end;
title;
%mend;
%macro anova(datain,varlist,group); 
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne );%let var=%scan(&varlist.,&varnumber.);
proc anova data=&datain; class &group; model &var = &group; run; 
%let varnumber=%eval(&varnumber.+1); %end;
%mend;
%macro cat(dat,cat);proc freq data=&dat; tables &cat/norow nocol nopercent;run;%mend;
%macro con(dat,con);proc means data=&dat n nmiss min p1 p25 p50 p75 p99 max; var &con; run;%mend;
/*get levels*/
%macro getlevels(dataset, var);
title "&dataset and &var";
proc sql; select count(distinct &var) as Levels,
count(*) as Nobs
from &dataset;
quit;
title;
%mend;

/*define per SD*/
%macro perSD(datain,varlist);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne ); %let var=%scan(&varlist.,&varnumber.);
/*standardization of exposure*/
proc means data=&datain; var &var; output out=stdfile; run;
proc sql; create table tmp as select &var as std from stdfile where _stat_="STD"; quit;
data &datain; if _n_=1 then set tmp (keep=std); set &datain; if std in (0,.) then &var._s=.; else &var._s=&var/std; run;
%let varnumber=%eval(&varnumber.+1);
%end;
%mend;
%macro perUnit(datain,varlist,units);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne ); %let var=%scan(&varlist.,&varnumber.);
/*standardization of exposure*/
data &datain; set &datain; &var._u=&var/&units; run;
%let varnumber=%eval(&varnumber.+1);
%end;
%mend;

/*assign median value for trend test*/
%macro median(datain,varlist,suffix=t);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne ); %let var=%scan(&varlist.,&varnumber.);
/*assign median value for each percentile*/
proc sort data=&datain; by &var.&suffix;run;
proc means data=&datain; var &var;class &var.&suffix; output out=medi median=median; run; 
data media; set medi;
if &var.&suffix=. then delete;
&var._m=median;
keep &var.&suffix &var._m;
rename &var.&suffix=group;
run;
data mlinear;set media; rename group=&var.&suffix;run;
data &datain; merge &datain mlinear; by &var.&suffix; run;
%let varnumber=%eval(&varnumber.+1);
%end;
%mend;

%macro checkvar (datain,dec,cov);
proc means n nmiss min p1 p25 median p75 p99 max mean std data=&datain maxdec=&dec; var &cov; run;
%mend;

libname datloc "/udd/nhlon/ukbssb/UKBdata";
options mautosource sasautos="/usr/local/channing/sasautos";


/*************************************************************************************************************************/
/*                                              Part I : Variables readin-proteomics                                     */
/*************************************************************************************************************************/

/*************************************************************************************************************************/
/*                                              Part I : Variables readin-metabolomics                                   */
/*************************************************************************************************************************/

/*************************************************************************************************************************/
/*                                              Part I : Variables readin-exposure                                       */
/*************************************************************************************************************************/
data hpyloriold; set datloc.hpylori;
if n_23039_0_0=. then  CagA=.; else if n_23039_0_0>400 then  CagA=1; else  CagA=0; /*definition I*/
if n_23040_0_0=. then  VacA=.; else if n_23040_0_0>100 then  VacA=1; else  VacA=0; /*definition II*/
if n_23041_0_0=. then   OMP=.; else if n_23041_0_0>170 then   OMP=1; else   OMP=0;
if n_23042_0_0=. then GroEL=.; else if n_23042_0_0> 80 then GroEL=1; else GroEL=0;
if n_23043_0_0=. then  Cata=.; else if n_23043_0_0>180 then  Cata=1; else  Cata=0;
if n_23044_0_0=. then  UreA=.; else if n_23044_0_0>130 then  UreA=1; else  UreA=0;

sumhpy1=sum(CagA,VacA,OMP,GroEL,Cata,UreA);
if sumhpy1=. then hpylori1=.;
 else if sumhpy1>=2 then hpylori1=1;
  else                   hpylori1=0;

sumhpy2=sum(VacA,OMP,GroEL,Cata,UreA);
if sumhpy2=. then hpylori2=.;
 else if sumhpy2>=2 then hpylori2=1;
  else                   hpylori2=0;

ind_QC=n_23049_0_0;
/*starttdate=s_23048_0_0;*/
run;

libname blood "/udd/nhlon/ukbssb/UKBdata";
data bloodday; set blood.blooddate; run;

%getlevels(hpyloriold, n_eid);
%getlevels(bloodday, n_eid);

proc sort data=hpyloriold; by n_eid; run;
proc sort data=bloodday;   by n_eid; run;

data hpylori; merge hpyloriold(in=a) bloodday(in=b); by n_eid; if a and b; 
starttdate = datepart(s_3166_0_0);
run;

proc freq data=hpylori;
 table VacA OMP GroEL Cata UreA sumhpy1 sumhpy2 hpylori1 hpylori2 ind_QC;
run;

proc means n nmiss min p25 p50 p75 max mean std data=hpylori;
 var n_23039_0_0 n_23039_1_0
     n_23040_0_0 n_23040_1_0
     n_23041_0_0 n_23041_1_0
     n_23042_0_0 n_23042_1_0
     n_23043_0_0 n_23043_1_0
     n_23044_0_0 n_23044_1_0;
run;

/*************************************************************************************************************************/
/*                                              Part I : Variables readin-outcome                                        */
/*************************************************************************************************************************/
data outcome; set datloc.outcomeset;
/*Define liver cancer from cancer registry*/
array ICD10 {*} s_40006_0_0 s_40006_1_0 s_40006_2_0 s_40006_3_0 s_40006_4_0 s_40006_5_0 s_40006_6_0 s_40006_7_0 
                s_40006_8_0 s_40006_9_0 s_40006_10_0 s_40006_11_0 s_40006_12_0 s_40006_13_0 s_40006_14_0 s_40006_15_0 
                s_40006_16_0 s_40006_17_0 s_40006_18_0 s_40006_19_0 s_40006_20_0 s_40006_21_0;
array cadate {*} s_40005_0_0 s_40005_1_0 s_40005_2_0 s_40005_3_0 s_40005_4_0 s_40005_5_0 s_40005_6_0 s_40005_7_0 
                 s_40005_8_0 s_40005_9_0 s_40005_10_0 s_40005_11_0 s_40005_12_0 s_40005_13_0 s_40005_14_0 s_40005_15_0 
                 s_40005_16_0 s_40005_17_0 s_40005_18_0 s_40005_19_0 s_40005_20_0 s_40005_21_0;

array liCdate {*} liCdate0-LiCdate21;
array HCCdate {*} HCCdate0-HCCdate21;
array ICCdate {*} ICCdate0-ICCdate21;

care_liver=0;care_HCC=0;care_ICC=0;
do i=1 to dim(ICD10);
if ICD10{i} in ("C220","C221","C223","C224","C227","C229") then do; care_liver=1; liCdate{i}=cadate{i};end; else do; liCdate{i}=.;end;
if ICD10{i} in ("C220")                                    then do;   care_HCC=1; HCCdate{i}=cadate{i};end; else do; HCCdate{i}=.;end;
if ICD10{i} in ("C221")                                    then do;   care_ICC=1; ICCdate{i}=cadate{i};end; else do; ICCdate{i}=.;end;
end;

caredate_lic=min(liCdate0, liCdate1, liCdate2, liCdate3, liCdate4, liCdate5, liCdate6, liCdate7, liCdate8, liCdate9, 
                 liCdate10, liCdate11, liCdate12, liCdate13, liCdate14, liCdate15, liCdate16, liCdate17, liCdate18, 
                 liCdate19, liCdate20, liCdate21);

caredate_hcc=min( HCCdate0,  HCCdate1, HCCdate2, HCCdate3, HCCdate4, HCCdate5, HCCdate6, HCCdate7, HCCdate8, HCCdate9, 
                  HCCdate10,  HCCdate11, HCCdate12, HCCdate13, HCCdate14, HCCdate15, HCCdate16, HCCdate17, HCCdate18, 
                  HCCdate19,  HCCdate20, HCCdate21);

caredate_icc=min( ICCdate0, ICCdate1, ICCdate2, ICCdate3, ICCdate4, ICCdate5, ICCdate6, ICCdate7, ICCdate8, ICCdate9, 
                  ICCdate10, ICCdate11, ICCdate12, ICCdate13, ICCdate14, ICCdate15, ICCdate16, ICCdate17, ICCdate18, 
                  ICCdate19, ICCdate20, ICCdate21);

/*Define liver disease from hospital data*/
array hspicd10 {*} s_41270_0_0 s_41270_0_1 s_41270_0_2 s_41270_0_3 s_41270_0_4 s_41270_0_5 s_41270_0_6 s_41270_0_7 s_41270_0_8 s_41270_0_9 s_41270_0_10 
                  s_41270_0_11 s_41270_0_12 s_41270_0_13 s_41270_0_14 s_41270_0_15 s_41270_0_16 s_41270_0_17 s_41270_0_18 s_41270_0_19 s_41270_0_20 
                  s_41270_0_21 s_41270_0_22 s_41270_0_23 s_41270_0_24 s_41270_0_25 s_41270_0_26 s_41270_0_27 s_41270_0_28 s_41270_0_29 s_41270_0_30 
                  s_41270_0_31 s_41270_0_32 s_41270_0_33 s_41270_0_34 s_41270_0_35 s_41270_0_36 s_41270_0_37 s_41270_0_38 s_41270_0_39 s_41270_0_40 
                  s_41270_0_41 s_41270_0_42 s_41270_0_43 s_41270_0_44 s_41270_0_45 s_41270_0_46 s_41270_0_47 s_41270_0_48 s_41270_0_49 s_41270_0_50 
                  s_41270_0_51 s_41270_0_52 s_41270_0_53 s_41270_0_54 s_41270_0_55 s_41270_0_56 s_41270_0_57 s_41270_0_58 s_41270_0_59 s_41270_0_60 
                  s_41270_0_61 s_41270_0_62 s_41270_0_63 s_41270_0_64 s_41270_0_65 s_41270_0_66 s_41270_0_67 s_41270_0_68 s_41270_0_69 s_41270_0_70 
                  s_41270_0_71 s_41270_0_72 s_41270_0_73 s_41270_0_74 s_41270_0_75 s_41270_0_76 s_41270_0_77 s_41270_0_78 s_41270_0_79 s_41270_0_80 
                  s_41270_0_81 s_41270_0_82 s_41270_0_83 s_41270_0_84 s_41270_0_85 s_41270_0_86 s_41270_0_87 s_41270_0_88 s_41270_0_89 s_41270_0_90 
                  s_41270_0_91 s_41270_0_92 s_41270_0_93 s_41270_0_94 s_41270_0_95 s_41270_0_96 s_41270_0_97 s_41270_0_98 s_41270_0_99 s_41270_0_100 
                  s_41270_0_101 s_41270_0_102 s_41270_0_103 s_41270_0_104 s_41270_0_105 s_41270_0_106 s_41270_0_107 s_41270_0_108 s_41270_0_109 s_41270_0_110 
                  s_41270_0_111 s_41270_0_112 s_41270_0_113 s_41270_0_114 s_41270_0_115 s_41270_0_116 s_41270_0_117 s_41270_0_118 s_41270_0_119 s_41270_0_120 
                  s_41270_0_121 s_41270_0_122 s_41270_0_123 s_41270_0_124 s_41270_0_125 s_41270_0_126 s_41270_0_127 s_41270_0_128 s_41270_0_129 s_41270_0_130 
                  s_41270_0_131 s_41270_0_132 s_41270_0_133 s_41270_0_134 s_41270_0_135 s_41270_0_136 s_41270_0_137 s_41270_0_138 s_41270_0_139 s_41270_0_140 
                  s_41270_0_141 s_41270_0_142 s_41270_0_143 s_41270_0_144 s_41270_0_145 s_41270_0_146 s_41270_0_147 s_41270_0_148 s_41270_0_149 s_41270_0_150 
                  s_41270_0_151 s_41270_0_152 s_41270_0_153 s_41270_0_154 s_41270_0_155 s_41270_0_156 s_41270_0_157 s_41270_0_158 s_41270_0_159 s_41270_0_160 
                  s_41270_0_161 s_41270_0_162 s_41270_0_163 s_41270_0_164 s_41270_0_165 s_41270_0_166 s_41270_0_167 s_41270_0_168 s_41270_0_169 s_41270_0_170 
                  s_41270_0_171 s_41270_0_172 s_41270_0_173 s_41270_0_174 s_41270_0_175 s_41270_0_176 s_41270_0_177 s_41270_0_178 s_41270_0_179 s_41270_0_180 
                  s_41270_0_181 s_41270_0_182 s_41270_0_183 s_41270_0_184 s_41270_0_185 s_41270_0_186 s_41270_0_187 s_41270_0_188 s_41270_0_189 s_41270_0_190 
                  s_41270_0_191 s_41270_0_192 s_41270_0_193 s_41270_0_194 s_41270_0_195 s_41270_0_196 s_41270_0_197 s_41270_0_198 s_41270_0_199 s_41270_0_200 
                  s_41270_0_201 s_41270_0_202 s_41270_0_203 s_41270_0_204 s_41270_0_205 s_41270_0_206 s_41270_0_207 s_41270_0_208 s_41270_0_209 s_41270_0_210 
                  s_41270_0_211 s_41270_0_212 s_41270_0_213 s_41270_0_214 s_41270_0_215 s_41270_0_216 s_41270_0_217 s_41270_0_218 s_41270_0_219 s_41270_0_220 
                  s_41270_0_221 s_41270_0_222 s_41270_0_223 s_41270_0_224 s_41270_0_225 s_41270_0_226 s_41270_0_227 s_41270_0_228 s_41270_0_229 s_41270_0_230 
                  s_41270_0_231 s_41270_0_232 s_41270_0_233 s_41270_0_234 s_41270_0_235 s_41270_0_236 s_41270_0_237 s_41270_0_238 s_41270_0_239 s_41270_0_240 
                  s_41270_0_241 s_41270_0_242;
array hspidt {*} s_41280_0_0 s_41280_0_1 s_41280_0_2 s_41280_0_3 s_41280_0_4 s_41280_0_5 s_41280_0_6 s_41280_0_7 s_41280_0_8 s_41280_0_9 s_41280_0_10 
                  s_41280_0_11 s_41280_0_12 s_41280_0_13 s_41280_0_14 s_41280_0_15 s_41280_0_16 s_41280_0_17 s_41280_0_18 s_41280_0_19 s_41280_0_20 
                  s_41280_0_21 s_41280_0_22 s_41280_0_23 s_41280_0_24 s_41280_0_25 s_41280_0_26 s_41280_0_27 s_41280_0_28 s_41280_0_29 s_41280_0_30 
                  s_41280_0_31 s_41280_0_32 s_41280_0_33 s_41280_0_34 s_41280_0_35 s_41280_0_36 s_41280_0_37 s_41280_0_38 s_41280_0_39 s_41280_0_40 
                  s_41280_0_41 s_41280_0_42 s_41280_0_43 s_41280_0_44 s_41280_0_45 s_41280_0_46 s_41280_0_47 s_41280_0_48 s_41280_0_49 s_41280_0_50 
                  s_41280_0_51 s_41280_0_52 s_41280_0_53 s_41280_0_54 s_41280_0_55 s_41280_0_56 s_41280_0_57 s_41280_0_58 s_41280_0_59 s_41280_0_60 
                  s_41280_0_61 s_41280_0_62 s_41280_0_63 s_41280_0_64 s_41280_0_65 s_41280_0_66 s_41280_0_67 s_41280_0_68 s_41280_0_69 s_41280_0_70 
                  s_41280_0_71 s_41280_0_72 s_41280_0_73 s_41280_0_74 s_41280_0_75 s_41280_0_76 s_41280_0_77 s_41280_0_78 s_41280_0_79 s_41280_0_80 
                  s_41280_0_81 s_41280_0_82 s_41280_0_83 s_41280_0_84 s_41280_0_85 s_41280_0_86 s_41280_0_87 s_41280_0_88 s_41280_0_89 s_41280_0_90 
                  s_41280_0_91 s_41280_0_92 s_41280_0_93 s_41280_0_94 s_41280_0_95 s_41280_0_96 s_41280_0_97 s_41280_0_98 s_41280_0_99 s_41280_0_100 
                  s_41280_0_101 s_41280_0_102 s_41280_0_103 s_41280_0_104 s_41280_0_105 s_41280_0_106 s_41280_0_107 s_41280_0_108 s_41280_0_109 s_41280_0_110 
                  s_41280_0_111 s_41280_0_112 s_41280_0_113 s_41280_0_114 s_41280_0_115 s_41280_0_116 s_41280_0_117 s_41280_0_118 s_41280_0_119 s_41280_0_120 
                  s_41280_0_121 s_41280_0_122 s_41280_0_123 s_41280_0_124 s_41280_0_125 s_41280_0_126 s_41280_0_127 s_41280_0_128 s_41280_0_129 s_41280_0_130 
                  s_41280_0_131 s_41280_0_132 s_41280_0_133 s_41280_0_134 s_41280_0_135 s_41280_0_136 s_41280_0_137 s_41280_0_138 s_41280_0_139 s_41280_0_140 
                  s_41280_0_141 s_41280_0_142 s_41280_0_143 s_41280_0_144 s_41280_0_145 s_41280_0_146 s_41280_0_147 s_41280_0_148 s_41280_0_149 s_41280_0_150 
                  s_41280_0_151 s_41280_0_152 s_41280_0_153 s_41280_0_154 s_41280_0_155 s_41280_0_156 s_41280_0_157 s_41280_0_158 s_41280_0_159 s_41280_0_160 
                  s_41280_0_161 s_41280_0_162 s_41280_0_163 s_41280_0_164 s_41280_0_165 s_41280_0_166 s_41280_0_167 s_41280_0_168 s_41280_0_169 s_41280_0_170 
                  s_41280_0_171 s_41280_0_172 s_41280_0_173 s_41280_0_174 s_41280_0_175 s_41280_0_176 s_41280_0_177 s_41280_0_178 s_41280_0_179 s_41280_0_180 
                  s_41280_0_181 s_41280_0_182 s_41280_0_183 s_41280_0_184 s_41280_0_185 s_41280_0_186 s_41280_0_187 s_41280_0_188 s_41280_0_189 s_41280_0_190 
                  s_41280_0_191 s_41280_0_192 s_41280_0_193 s_41280_0_194 s_41280_0_195 s_41280_0_196 s_41280_0_197 s_41280_0_198 s_41280_0_199 s_41280_0_200 
                  s_41280_0_201 s_41280_0_202 s_41280_0_203 s_41280_0_204 s_41280_0_205 s_41280_0_206 s_41280_0_207 s_41280_0_208 s_41280_0_209 s_41280_0_210 
                  s_41280_0_211 s_41280_0_212 s_41280_0_213 s_41280_0_214 s_41280_0_215 s_41280_0_216 s_41280_0_217 s_41280_0_218 s_41280_0_219 s_41280_0_220 
                  s_41280_0_221 s_41280_0_222 s_41280_0_223 s_41280_0_224 s_41280_0_225 s_41280_0_226 s_41280_0_227 s_41280_0_228 s_41280_0_229 s_41280_0_230 
                  s_41280_0_231 s_41280_0_232 s_41280_0_233 s_41280_0_234 s_41280_0_235 s_41280_0_236 s_41280_0_237 s_41280_0_238 s_41280_0_239 s_41280_0_240 
                  s_41280_0_241 s_41280_0_242;
array naflddt     {*} naflddt0-naflddt242;
array cirrhosisdt {*} cirrhosisdt0-cirrhosisdt242;
array crhs_cpdt   {*} crhs_cpdt0-crhs_cpdt242;
array crhs_dcdt   {*} crhs_dcdt0-crhs_dcdt242;
array slddt       {*} slddt0-slddt242; 
array lic_hsdt    {*} lic_hsdt0-lic_hsdt242;
array hcc_hsdt    {*} hcc_hsdt0-hcc_hsdt242;
array icc_hsdt    {*} icc_hsdt0-icc_hsdt242;

nafld=0; cirrhosis=0; sld=0; crhs_cp=0; crhs_dc=0; lic_hs=0; HCC_hs=0; ICC_hs=0;

do i=1 to dim(hspicd10);
if hspicd10{i} in ("K758","K760")                                    then do; nafld=1;     naflddt{i}=hspidt{i};    end;else do; naflddt{i}=.;    end;
if hspicd10{i} in ("K740","K741","K742","K743","K744","K745","K746") then do; cirrhosis=1; cirrhosisdt{i}=hspidt{i};end;else do; cirrhosisdt{i}=.;end;
if hspicd10{i} in ("K746","I859","I982","I864")                      then do; crhs_cp=1;   crhs_cpdt{i}=hspidt{i};  end;else do; crhs_cpdt{i}=.;  end;
if hspicd10{i} in ("I850","I983","R18","K767","K766")                then do; crhs_dc=1;   crhs_dcdt{i}=hspidt{i};  end;else do; crhs_dcdt{i}=.;  end;
if hspicd10{i} in ("C220","I850","I859","K703","K704","K721","K729",
                   "K741","K742","K746","K766","K767","Z944")        then do; sld=1;       slddt{i}=hspidt{i};      end;else do; slddt{i}=.;      end;
if hspicd10{i} in ("C220","C221","C223","C224","C227","C229")        then do; lic_hs=1;    lic_hsdt{i}=hspidt{i};   end;else do; lic_hsdt{i}=.;   end;
if hspicd10{i} in ("C220")                                           then do; hcc_hs=1;    hcc_hsdt{i}=hspidt{i};   end;else do; hcc_hsdt{i}=.;   end;
if hspicd10{i} in ("C221")                                           then do; icc_hs=1;    icc_hsdt{i}=hspidt{i};   end;else do; icc_hsdt{i}=.;   end;
/*compensated cirrhosis and decompensated cirrhosis based on 10.1002/hep.31726, Hagstrom, 2021, Hepatology*/
end;

naflddate=min(naflddt0, naflddt1, naflddt2, naflddt3, naflddt4, naflddt5, naflddt6, naflddt7, naflddt8, naflddt9, naflddt10, naflddt11, naflddt12, naflddt13, naflddt14, naflddt15, naflddt16, naflddt17, naflddt18, naflddt19, naflddt20, naflddt21, naflddt22, naflddt23, naflddt24, naflddt25, naflddt26, naflddt27, naflddt28, naflddt29, naflddt30, naflddt31, naflddt32, naflddt33, naflddt34, naflddt35, naflddt36, naflddt37, naflddt38, naflddt39, naflddt40, naflddt41, naflddt42, naflddt43, naflddt44, naflddt45, naflddt46, naflddt47, naflddt48, naflddt49, naflddt50, naflddt51, naflddt52, naflddt53, naflddt54, naflddt55, naflddt56, naflddt57, naflddt58, naflddt59, naflddt60, naflddt61, naflddt62, naflddt63, naflddt64, naflddt65, naflddt66, naflddt67, naflddt68, naflddt69, naflddt70, naflddt71, naflddt72, naflddt73, naflddt74, naflddt75, naflddt76, naflddt77, naflddt78, naflddt79, naflddt80, naflddt81, naflddt82, naflddt83, naflddt84, naflddt85, naflddt86, naflddt87, naflddt88, naflddt89, naflddt90, naflddt91, naflddt92, naflddt93, naflddt94, naflddt95, naflddt96, naflddt97, naflddt98, naflddt99, naflddt100, naflddt101, naflddt102, naflddt103, naflddt104, naflddt105, naflddt106, naflddt107, naflddt108, naflddt109, naflddt110, naflddt111, naflddt112, naflddt113, naflddt114, naflddt115, naflddt116, naflddt117, naflddt118, naflddt119, naflddt120, naflddt121, naflddt122, naflddt123, naflddt124, naflddt125, naflddt126, naflddt127, naflddt128, naflddt129, naflddt130, naflddt131, naflddt132, naflddt133, naflddt134, naflddt135, naflddt136, naflddt137, naflddt138, naflddt139, naflddt140, naflddt141, naflddt142, naflddt143, naflddt144, naflddt145, naflddt146, naflddt147, naflddt148, naflddt149, naflddt150, naflddt151, naflddt152, naflddt153, naflddt154, naflddt155, naflddt156, naflddt157, naflddt158, naflddt159, naflddt160, naflddt161, naflddt162, naflddt163, naflddt164, naflddt165, naflddt166, naflddt167, naflddt168, naflddt169, naflddt170, naflddt171, naflddt172, naflddt173, naflddt174, naflddt175, naflddt176, naflddt177, naflddt178, naflddt179, naflddt180, naflddt181, naflddt182, naflddt183, naflddt184, naflddt185, naflddt186, naflddt187, naflddt188, naflddt189, naflddt190, naflddt191, naflddt192, naflddt193, naflddt194, naflddt195, naflddt196, naflddt197, naflddt198, naflddt199, naflddt200, naflddt201, naflddt202, naflddt203, naflddt204, naflddt205, naflddt206, naflddt207, naflddt208, naflddt209, naflddt210, naflddt211, naflddt212, naflddt213, naflddt214, naflddt215, naflddt216, naflddt217, naflddt218, naflddt219, naflddt220, naflddt221, naflddt222, naflddt223, naflddt224, naflddt225, naflddt226, naflddt227, naflddt228, naflddt229, naflddt230, naflddt231, naflddt232, naflddt233, naflddt234, naflddt235, naflddt236, naflddt237, naflddt238, naflddt239, naflddt240, naflddt241, naflddt242);
cirrhosisdate=min(cirrhosisdt0, cirrhosisdt1, cirrhosisdt2, cirrhosisdt3, cirrhosisdt4, cirrhosisdt5, cirrhosisdt6, cirrhosisdt7, cirrhosisdt8, cirrhosisdt9, cirrhosisdt10, cirrhosisdt11, cirrhosisdt12, cirrhosisdt13, cirrhosisdt14, cirrhosisdt15, cirrhosisdt16, cirrhosisdt17, cirrhosisdt18, cirrhosisdt19, cirrhosisdt20, cirrhosisdt21, cirrhosisdt22, cirrhosisdt23, cirrhosisdt24, cirrhosisdt25, cirrhosisdt26, cirrhosisdt27, cirrhosisdt28, cirrhosisdt29, cirrhosisdt30, cirrhosisdt31, cirrhosisdt32, cirrhosisdt33, cirrhosisdt34, cirrhosisdt35, cirrhosisdt36, cirrhosisdt37, cirrhosisdt38, cirrhosisdt39, cirrhosisdt40, cirrhosisdt41, cirrhosisdt42, cirrhosisdt43, cirrhosisdt44, cirrhosisdt45, cirrhosisdt46, cirrhosisdt47, cirrhosisdt48, cirrhosisdt49, cirrhosisdt50, cirrhosisdt51, cirrhosisdt52, cirrhosisdt53, cirrhosisdt54, cirrhosisdt55, cirrhosisdt56, cirrhosisdt57, cirrhosisdt58, cirrhosisdt59, cirrhosisdt60, cirrhosisdt61, cirrhosisdt62, cirrhosisdt63, cirrhosisdt64, cirrhosisdt65, cirrhosisdt66, cirrhosisdt67, cirrhosisdt68, cirrhosisdt69, cirrhosisdt70, cirrhosisdt71, cirrhosisdt72, cirrhosisdt73, cirrhosisdt74, cirrhosisdt75, cirrhosisdt76, cirrhosisdt77, cirrhosisdt78, cirrhosisdt79, cirrhosisdt80, cirrhosisdt81, cirrhosisdt82, cirrhosisdt83, cirrhosisdt84, cirrhosisdt85, cirrhosisdt86, cirrhosisdt87, cirrhosisdt88, cirrhosisdt89, cirrhosisdt90, cirrhosisdt91, cirrhosisdt92, cirrhosisdt93, cirrhosisdt94, cirrhosisdt95, cirrhosisdt96, cirrhosisdt97, cirrhosisdt98, cirrhosisdt99, cirrhosisdt100, cirrhosisdt101, cirrhosisdt102, cirrhosisdt103, cirrhosisdt104, cirrhosisdt105, cirrhosisdt106, cirrhosisdt107, cirrhosisdt108, cirrhosisdt109, cirrhosisdt110, cirrhosisdt111, cirrhosisdt112, cirrhosisdt113, cirrhosisdt114, cirrhosisdt115, cirrhosisdt116, cirrhosisdt117, cirrhosisdt118, cirrhosisdt119, cirrhosisdt120, cirrhosisdt121, cirrhosisdt122, cirrhosisdt123, cirrhosisdt124, cirrhosisdt125, cirrhosisdt126, cirrhosisdt127, cirrhosisdt128, cirrhosisdt129, cirrhosisdt130, cirrhosisdt131, cirrhosisdt132, cirrhosisdt133, cirrhosisdt134, cirrhosisdt135, cirrhosisdt136, cirrhosisdt137, cirrhosisdt138, cirrhosisdt139, cirrhosisdt140, cirrhosisdt141, cirrhosisdt142, cirrhosisdt143, cirrhosisdt144, cirrhosisdt145, cirrhosisdt146, cirrhosisdt147, cirrhosisdt148, cirrhosisdt149, cirrhosisdt150, cirrhosisdt151, cirrhosisdt152, cirrhosisdt153, cirrhosisdt154, cirrhosisdt155, cirrhosisdt156, cirrhosisdt157, cirrhosisdt158, cirrhosisdt159, cirrhosisdt160, cirrhosisdt161, cirrhosisdt162, cirrhosisdt163, cirrhosisdt164, cirrhosisdt165, cirrhosisdt166, cirrhosisdt167, cirrhosisdt168, cirrhosisdt169, cirrhosisdt170, cirrhosisdt171, cirrhosisdt172, cirrhosisdt173, cirrhosisdt174, cirrhosisdt175, cirrhosisdt176, cirrhosisdt177, cirrhosisdt178, cirrhosisdt179, cirrhosisdt180, cirrhosisdt181, cirrhosisdt182, cirrhosisdt183, cirrhosisdt184, cirrhosisdt185, cirrhosisdt186, cirrhosisdt187, cirrhosisdt188, cirrhosisdt189, cirrhosisdt190, cirrhosisdt191, cirrhosisdt192, cirrhosisdt193, cirrhosisdt194, cirrhosisdt195, cirrhosisdt196, cirrhosisdt197, cirrhosisdt198, cirrhosisdt199, cirrhosisdt200, cirrhosisdt201, cirrhosisdt202, cirrhosisdt203, cirrhosisdt204, cirrhosisdt205, cirrhosisdt206, cirrhosisdt207, cirrhosisdt208, cirrhosisdt209, cirrhosisdt210, cirrhosisdt211, cirrhosisdt212, cirrhosisdt213, cirrhosisdt214, cirrhosisdt215, cirrhosisdt216, cirrhosisdt217, cirrhosisdt218, cirrhosisdt219, cirrhosisdt220, cirrhosisdt221, cirrhosisdt222, cirrhosisdt223, cirrhosisdt224, cirrhosisdt225, cirrhosisdt226, cirrhosisdt227, cirrhosisdt228, cirrhosisdt229, cirrhosisdt230, cirrhosisdt231, cirrhosisdt232, cirrhosisdt233, cirrhosisdt234, cirrhosisdt235, cirrhosisdt236, cirrhosisdt237, cirrhosisdt238, cirrhosisdt239, cirrhosisdt240, cirrhosisdt241, cirrhosisdt242);
crhs_cpdate=min(crhs_cpdt0, crhs_cpdt1, crhs_cpdt2, crhs_cpdt3, crhs_cpdt4, crhs_cpdt5, crhs_cpdt6, crhs_cpdt7, crhs_cpdt8, crhs_cpdt9, crhs_cpdt10, crhs_cpdt11, crhs_cpdt12, crhs_cpdt13, crhs_cpdt14, crhs_cpdt15, crhs_cpdt16, crhs_cpdt17, crhs_cpdt18, crhs_cpdt19, crhs_cpdt20, crhs_cpdt21, crhs_cpdt22, crhs_cpdt23, crhs_cpdt24, crhs_cpdt25, crhs_cpdt26, crhs_cpdt27, crhs_cpdt28, crhs_cpdt29, crhs_cpdt30, crhs_cpdt31, crhs_cpdt32, crhs_cpdt33, crhs_cpdt34, crhs_cpdt35, crhs_cpdt36, crhs_cpdt37, crhs_cpdt38, crhs_cpdt39, crhs_cpdt40, crhs_cpdt41, crhs_cpdt42, crhs_cpdt43, crhs_cpdt44, crhs_cpdt45, crhs_cpdt46, crhs_cpdt47, crhs_cpdt48, crhs_cpdt49, crhs_cpdt50, crhs_cpdt51, crhs_cpdt52, crhs_cpdt53, crhs_cpdt54, crhs_cpdt55, crhs_cpdt56, crhs_cpdt57, crhs_cpdt58, crhs_cpdt59, crhs_cpdt60, crhs_cpdt61, crhs_cpdt62, crhs_cpdt63, crhs_cpdt64, crhs_cpdt65, crhs_cpdt66, crhs_cpdt67, crhs_cpdt68, crhs_cpdt69, crhs_cpdt70, crhs_cpdt71, crhs_cpdt72, crhs_cpdt73, crhs_cpdt74, crhs_cpdt75, crhs_cpdt76, crhs_cpdt77, crhs_cpdt78, crhs_cpdt79, crhs_cpdt80, crhs_cpdt81, crhs_cpdt82, crhs_cpdt83, crhs_cpdt84, crhs_cpdt85, crhs_cpdt86, crhs_cpdt87, crhs_cpdt88, crhs_cpdt89, crhs_cpdt90, crhs_cpdt91, crhs_cpdt92, crhs_cpdt93, crhs_cpdt94, crhs_cpdt95, crhs_cpdt96, crhs_cpdt97, crhs_cpdt98, crhs_cpdt99, crhs_cpdt100, crhs_cpdt101, crhs_cpdt102, crhs_cpdt103, crhs_cpdt104, crhs_cpdt105, crhs_cpdt106, crhs_cpdt107, crhs_cpdt108, crhs_cpdt109, crhs_cpdt110, crhs_cpdt111, crhs_cpdt112, crhs_cpdt113, crhs_cpdt114, crhs_cpdt115, crhs_cpdt116, crhs_cpdt117, crhs_cpdt118, crhs_cpdt119, crhs_cpdt120, crhs_cpdt121, crhs_cpdt122, crhs_cpdt123, crhs_cpdt124, crhs_cpdt125, crhs_cpdt126, crhs_cpdt127, crhs_cpdt128, crhs_cpdt129, crhs_cpdt130, crhs_cpdt131, crhs_cpdt132, crhs_cpdt133, crhs_cpdt134, crhs_cpdt135, crhs_cpdt136, crhs_cpdt137, crhs_cpdt138, crhs_cpdt139, crhs_cpdt140, crhs_cpdt141, crhs_cpdt142, crhs_cpdt143, crhs_cpdt144, crhs_cpdt145, crhs_cpdt146, crhs_cpdt147, crhs_cpdt148, crhs_cpdt149, crhs_cpdt150, crhs_cpdt151, crhs_cpdt152, crhs_cpdt153, crhs_cpdt154, crhs_cpdt155, crhs_cpdt156, crhs_cpdt157, crhs_cpdt158, crhs_cpdt159, crhs_cpdt160, crhs_cpdt161, crhs_cpdt162, crhs_cpdt163, crhs_cpdt164, crhs_cpdt165, crhs_cpdt166, crhs_cpdt167, crhs_cpdt168, crhs_cpdt169, crhs_cpdt170, crhs_cpdt171, crhs_cpdt172, crhs_cpdt173, crhs_cpdt174, crhs_cpdt175, crhs_cpdt176, crhs_cpdt177, crhs_cpdt178, crhs_cpdt179, crhs_cpdt180, crhs_cpdt181, crhs_cpdt182, crhs_cpdt183, crhs_cpdt184, crhs_cpdt185, crhs_cpdt186, crhs_cpdt187, crhs_cpdt188, crhs_cpdt189, crhs_cpdt190, crhs_cpdt191, crhs_cpdt192, crhs_cpdt193, crhs_cpdt194, crhs_cpdt195, crhs_cpdt196, crhs_cpdt197, crhs_cpdt198, crhs_cpdt199, crhs_cpdt200, crhs_cpdt201, crhs_cpdt202, crhs_cpdt203, crhs_cpdt204, crhs_cpdt205, crhs_cpdt206, crhs_cpdt207, crhs_cpdt208, crhs_cpdt209, crhs_cpdt210, crhs_cpdt211, crhs_cpdt212, crhs_cpdt213, crhs_cpdt214, crhs_cpdt215, crhs_cpdt216, crhs_cpdt217, crhs_cpdt218, crhs_cpdt219, crhs_cpdt220, crhs_cpdt221, crhs_cpdt222, crhs_cpdt223, crhs_cpdt224, crhs_cpdt225, crhs_cpdt226, crhs_cpdt227, crhs_cpdt228, crhs_cpdt229, crhs_cpdt230, crhs_cpdt231, crhs_cpdt232, crhs_cpdt233, crhs_cpdt234, crhs_cpdt235, crhs_cpdt236, crhs_cpdt237, crhs_cpdt238, crhs_cpdt239, crhs_cpdt240, crhs_cpdt241, crhs_cpdt242);
crhs_dcdate=min(crhs_dcdt0, crhs_dcdt1, crhs_dcdt2, crhs_dcdt3, crhs_dcdt4, crhs_dcdt5, crhs_dcdt6, crhs_dcdt7, crhs_dcdt8, crhs_dcdt9, crhs_dcdt10, crhs_dcdt11, crhs_dcdt12, crhs_dcdt13, crhs_dcdt14, crhs_dcdt15, crhs_dcdt16, crhs_dcdt17, crhs_dcdt18, crhs_dcdt19, crhs_dcdt20, crhs_dcdt21, crhs_dcdt22, crhs_dcdt23, crhs_dcdt24, crhs_dcdt25, crhs_dcdt26, crhs_dcdt27, crhs_dcdt28, crhs_dcdt29, crhs_dcdt30, crhs_dcdt31, crhs_dcdt32, crhs_dcdt33, crhs_dcdt34, crhs_dcdt35, crhs_dcdt36, crhs_dcdt37, crhs_dcdt38, crhs_dcdt39, crhs_dcdt40, crhs_dcdt41, crhs_dcdt42, crhs_dcdt43, crhs_dcdt44, crhs_dcdt45, crhs_dcdt46, crhs_dcdt47, crhs_dcdt48, crhs_dcdt49, crhs_dcdt50, crhs_dcdt51, crhs_dcdt52, crhs_dcdt53, crhs_dcdt54, crhs_dcdt55, crhs_dcdt56, crhs_dcdt57, crhs_dcdt58, crhs_dcdt59, crhs_dcdt60, crhs_dcdt61, crhs_dcdt62, crhs_dcdt63, crhs_dcdt64, crhs_dcdt65, crhs_dcdt66, crhs_dcdt67, crhs_dcdt68, crhs_dcdt69, crhs_dcdt70, crhs_dcdt71, crhs_dcdt72, crhs_dcdt73, crhs_dcdt74, crhs_dcdt75, crhs_dcdt76, crhs_dcdt77, crhs_dcdt78, crhs_dcdt79, crhs_dcdt80, crhs_dcdt81, crhs_dcdt82, crhs_dcdt83, crhs_dcdt84, crhs_dcdt85, crhs_dcdt86, crhs_dcdt87, crhs_dcdt88, crhs_dcdt89, crhs_dcdt90, crhs_dcdt91, crhs_dcdt92, crhs_dcdt93, crhs_dcdt94, crhs_dcdt95, crhs_dcdt96, crhs_dcdt97, crhs_dcdt98, crhs_dcdt99, crhs_dcdt100, crhs_dcdt101, crhs_dcdt102, crhs_dcdt103, crhs_dcdt104, crhs_dcdt105, crhs_dcdt106, crhs_dcdt107, crhs_dcdt108, crhs_dcdt109, crhs_dcdt110, crhs_dcdt111, crhs_dcdt112, crhs_dcdt113, crhs_dcdt114, crhs_dcdt115, crhs_dcdt116, crhs_dcdt117, crhs_dcdt118, crhs_dcdt119, crhs_dcdt120, crhs_dcdt121, crhs_dcdt122, crhs_dcdt123, crhs_dcdt124, crhs_dcdt125, crhs_dcdt126, crhs_dcdt127, crhs_dcdt128, crhs_dcdt129, crhs_dcdt130, crhs_dcdt131, crhs_dcdt132, crhs_dcdt133, crhs_dcdt134, crhs_dcdt135, crhs_dcdt136, crhs_dcdt137, crhs_dcdt138, crhs_dcdt139, crhs_dcdt140, crhs_dcdt141, crhs_dcdt142, crhs_dcdt143, crhs_dcdt144, crhs_dcdt145, crhs_dcdt146, crhs_dcdt147, crhs_dcdt148, crhs_dcdt149, crhs_dcdt150, crhs_dcdt151, crhs_dcdt152, crhs_dcdt153, crhs_dcdt154, crhs_dcdt155, crhs_dcdt156, crhs_dcdt157, crhs_dcdt158, crhs_dcdt159, crhs_dcdt160, crhs_dcdt161, crhs_dcdt162, crhs_dcdt163, crhs_dcdt164, crhs_dcdt165, crhs_dcdt166, crhs_dcdt167, crhs_dcdt168, crhs_dcdt169, crhs_dcdt170, crhs_dcdt171, crhs_dcdt172, crhs_dcdt173, crhs_dcdt174, crhs_dcdt175, crhs_dcdt176, crhs_dcdt177, crhs_dcdt178, crhs_dcdt179, crhs_dcdt180, crhs_dcdt181, crhs_dcdt182, crhs_dcdt183, crhs_dcdt184, crhs_dcdt185, crhs_dcdt186, crhs_dcdt187, crhs_dcdt188, crhs_dcdt189, crhs_dcdt190, crhs_dcdt191, crhs_dcdt192, crhs_dcdt193, crhs_dcdt194, crhs_dcdt195, crhs_dcdt196, crhs_dcdt197, crhs_dcdt198, crhs_dcdt199, crhs_dcdt200, crhs_dcdt201, crhs_dcdt202, crhs_dcdt203, crhs_dcdt204, crhs_dcdt205, crhs_dcdt206, crhs_dcdt207, crhs_dcdt208, crhs_dcdt209, crhs_dcdt210, crhs_dcdt211, crhs_dcdt212, crhs_dcdt213, crhs_dcdt214, crhs_dcdt215, crhs_dcdt216, crhs_dcdt217, crhs_dcdt218, crhs_dcdt219, crhs_dcdt220, crhs_dcdt221, crhs_dcdt222, crhs_dcdt223, crhs_dcdt224, crhs_dcdt225, crhs_dcdt226, crhs_dcdt227, crhs_dcdt228, crhs_dcdt229, crhs_dcdt230, crhs_dcdt231, crhs_dcdt232, crhs_dcdt233, crhs_dcdt234, crhs_dcdt235, crhs_dcdt236, crhs_dcdt237, crhs_dcdt238, crhs_dcdt239, crhs_dcdt240, crhs_dcdt241, crhs_dcdt242);
slddate=min(slddt0, slddt1, slddt2, slddt3, slddt4, slddt5, slddt6, slddt7, slddt8, slddt9, slddt10, slddt11, slddt12, slddt13, slddt14, slddt15, slddt16, slddt17, slddt18, slddt19, slddt20, slddt21, slddt22, slddt23, slddt24, slddt25, slddt26, slddt27, slddt28, slddt29, slddt30, slddt31, slddt32, slddt33, slddt34, slddt35, slddt36, slddt37, slddt38, slddt39, slddt40, slddt41, slddt42, slddt43, slddt44, slddt45, slddt46, slddt47, slddt48, slddt49, slddt50, slddt51, slddt52, slddt53, slddt54, slddt55, slddt56, slddt57, slddt58, slddt59, slddt60, slddt61, slddt62, slddt63, slddt64, slddt65, slddt66, slddt67, slddt68, slddt69, slddt70, slddt71, slddt72, slddt73, slddt74, slddt75, slddt76, slddt77, slddt78, slddt79, slddt80, slddt81, slddt82, slddt83, slddt84, slddt85, slddt86, slddt87, slddt88, slddt89, slddt90, slddt91, slddt92, slddt93, slddt94, slddt95, slddt96, slddt97, slddt98, slddt99, slddt100, slddt101, slddt102, slddt103, slddt104, slddt105, slddt106, slddt107, slddt108, slddt109, slddt110, slddt111, slddt112, slddt113, slddt114, slddt115, slddt116, slddt117, slddt118, slddt119, slddt120, slddt121, slddt122, slddt123, slddt124, slddt125, slddt126, slddt127, slddt128, slddt129, slddt130, slddt131, slddt132, slddt133, slddt134, slddt135, slddt136, slddt137, slddt138, slddt139, slddt140, slddt141, slddt142, slddt143, slddt144, slddt145, slddt146, slddt147, slddt148, slddt149, slddt150, slddt151, slddt152, slddt153, slddt154, slddt155, slddt156, slddt157, slddt158, slddt159, slddt160, slddt161, slddt162, slddt163, slddt164, slddt165, slddt166, slddt167, slddt168, slddt169, slddt170, slddt171, slddt172, slddt173, slddt174, slddt175, slddt176, slddt177, slddt178, slddt179, slddt180, slddt181, slddt182, slddt183, slddt184, slddt185, slddt186, slddt187, slddt188, slddt189, slddt190, slddt191, slddt192, slddt193, slddt194, slddt195, slddt196, slddt197, slddt198, slddt199, slddt200, slddt201, slddt202, slddt203, slddt204, slddt205, slddt206, slddt207, slddt208, slddt209, slddt210, slddt211, slddt212, slddt213, slddt214, slddt215, slddt216, slddt217, slddt218, slddt219, slddt220, slddt221, slddt222, slddt223, slddt224, slddt225, slddt226, slddt227, slddt228, slddt229, slddt230, slddt231, slddt232, slddt233, slddt234, slddt235, slddt236, slddt237, slddt238, slddt239, slddt240, slddt241, slddt242);
lic_hsdate=min(lic_hsdt0, lic_hsdt1, lic_hsdt2, lic_hsdt3, lic_hsdt4, lic_hsdt5, lic_hsdt6, lic_hsdt7, lic_hsdt8, lic_hsdt9, lic_hsdt10, lic_hsdt11, lic_hsdt12, lic_hsdt13, lic_hsdt14, lic_hsdt15, lic_hsdt16, lic_hsdt17, lic_hsdt18, lic_hsdt19, lic_hsdt20, lic_hsdt21, lic_hsdt22, lic_hsdt23, lic_hsdt24, lic_hsdt25, lic_hsdt26, lic_hsdt27, lic_hsdt28, lic_hsdt29, lic_hsdt30, lic_hsdt31, lic_hsdt32, lic_hsdt33, lic_hsdt34, lic_hsdt35, lic_hsdt36, lic_hsdt37, lic_hsdt38, lic_hsdt39, lic_hsdt40, lic_hsdt41, lic_hsdt42, lic_hsdt43, lic_hsdt44, lic_hsdt45, lic_hsdt46, lic_hsdt47, lic_hsdt48, lic_hsdt49, lic_hsdt50, lic_hsdt51, lic_hsdt52, lic_hsdt53, lic_hsdt54, lic_hsdt55, lic_hsdt56, lic_hsdt57, lic_hsdt58, lic_hsdt59, lic_hsdt60, lic_hsdt61, lic_hsdt62, lic_hsdt63, lic_hsdt64, lic_hsdt65, lic_hsdt66, lic_hsdt67, lic_hsdt68, lic_hsdt69, lic_hsdt70, lic_hsdt71, lic_hsdt72, lic_hsdt73, lic_hsdt74, lic_hsdt75, lic_hsdt76, lic_hsdt77, lic_hsdt78, lic_hsdt79, lic_hsdt80, lic_hsdt81, lic_hsdt82, lic_hsdt83, lic_hsdt84, lic_hsdt85, lic_hsdt86, lic_hsdt87, lic_hsdt88, lic_hsdt89, lic_hsdt90, lic_hsdt91, lic_hsdt92, lic_hsdt93, lic_hsdt94, lic_hsdt95, lic_hsdt96, lic_hsdt97, lic_hsdt98, lic_hsdt99, lic_hsdt100, lic_hsdt101, lic_hsdt102, lic_hsdt103, lic_hsdt104, lic_hsdt105, lic_hsdt106, lic_hsdt107, lic_hsdt108, lic_hsdt109, lic_hsdt110, lic_hsdt111, lic_hsdt112, lic_hsdt113, lic_hsdt114, lic_hsdt115, lic_hsdt116, lic_hsdt117, lic_hsdt118, lic_hsdt119, lic_hsdt120, lic_hsdt121, lic_hsdt122, lic_hsdt123, lic_hsdt124, lic_hsdt125, lic_hsdt126, lic_hsdt127, lic_hsdt128, lic_hsdt129, lic_hsdt130, lic_hsdt131, lic_hsdt132, lic_hsdt133, lic_hsdt134, lic_hsdt135, lic_hsdt136, lic_hsdt137, lic_hsdt138, lic_hsdt139, lic_hsdt140, lic_hsdt141, lic_hsdt142, lic_hsdt143, lic_hsdt144, lic_hsdt145, lic_hsdt146, lic_hsdt147, lic_hsdt148, lic_hsdt149, lic_hsdt150, lic_hsdt151, lic_hsdt152, lic_hsdt153, lic_hsdt154, lic_hsdt155, lic_hsdt156, lic_hsdt157, lic_hsdt158, lic_hsdt159, lic_hsdt160, lic_hsdt161, lic_hsdt162, lic_hsdt163, lic_hsdt164, lic_hsdt165, lic_hsdt166, lic_hsdt167, lic_hsdt168, lic_hsdt169, lic_hsdt170, lic_hsdt171, lic_hsdt172, lic_hsdt173, lic_hsdt174, lic_hsdt175, lic_hsdt176, lic_hsdt177, lic_hsdt178, lic_hsdt179, lic_hsdt180, lic_hsdt181, lic_hsdt182, lic_hsdt183, lic_hsdt184, lic_hsdt185, lic_hsdt186, lic_hsdt187, lic_hsdt188, lic_hsdt189, lic_hsdt190, lic_hsdt191, lic_hsdt192, lic_hsdt193, lic_hsdt194, lic_hsdt195, lic_hsdt196, lic_hsdt197, lic_hsdt198, lic_hsdt199, lic_hsdt200, lic_hsdt201, lic_hsdt202, lic_hsdt203, lic_hsdt204, lic_hsdt205, lic_hsdt206, lic_hsdt207, lic_hsdt208, lic_hsdt209, lic_hsdt210, lic_hsdt211, lic_hsdt212, lic_hsdt213, lic_hsdt214, lic_hsdt215, lic_hsdt216, lic_hsdt217, lic_hsdt218, lic_hsdt219, lic_hsdt220, lic_hsdt221, lic_hsdt222, lic_hsdt223, lic_hsdt224, lic_hsdt225, lic_hsdt226, lic_hsdt227, lic_hsdt228, lic_hsdt229, lic_hsdt230, lic_hsdt231, lic_hsdt232, lic_hsdt233, lic_hsdt234, lic_hsdt235, lic_hsdt236, lic_hsdt237, lic_hsdt238, lic_hsdt239, lic_hsdt240, lic_hsdt241, lic_hsdt242);
hcc_hsdate=min(hcc_hsdt0, hcc_hsdt1, hcc_hsdt2, hcc_hsdt3, hcc_hsdt4, hcc_hsdt5, hcc_hsdt6, hcc_hsdt7, hcc_hsdt8, hcc_hsdt9, hcc_hsdt10, hcc_hsdt11, hcc_hsdt12, hcc_hsdt13, hcc_hsdt14, hcc_hsdt15, hcc_hsdt16, hcc_hsdt17, hcc_hsdt18, hcc_hsdt19, hcc_hsdt20, hcc_hsdt21, hcc_hsdt22, hcc_hsdt23, hcc_hsdt24, hcc_hsdt25, hcc_hsdt26, hcc_hsdt27, hcc_hsdt28, hcc_hsdt29, hcc_hsdt30, hcc_hsdt31, hcc_hsdt32, hcc_hsdt33, hcc_hsdt34, hcc_hsdt35, hcc_hsdt36, hcc_hsdt37, hcc_hsdt38, hcc_hsdt39, hcc_hsdt40, hcc_hsdt41, hcc_hsdt42, hcc_hsdt43, hcc_hsdt44, hcc_hsdt45, hcc_hsdt46, hcc_hsdt47, hcc_hsdt48, hcc_hsdt49, hcc_hsdt50, hcc_hsdt51, hcc_hsdt52, hcc_hsdt53, hcc_hsdt54, hcc_hsdt55, hcc_hsdt56, hcc_hsdt57, hcc_hsdt58, hcc_hsdt59, hcc_hsdt60, hcc_hsdt61, hcc_hsdt62, hcc_hsdt63, hcc_hsdt64, hcc_hsdt65, hcc_hsdt66, hcc_hsdt67, hcc_hsdt68, hcc_hsdt69, hcc_hsdt70, hcc_hsdt71, hcc_hsdt72, hcc_hsdt73, hcc_hsdt74, hcc_hsdt75, hcc_hsdt76, hcc_hsdt77, hcc_hsdt78, hcc_hsdt79, hcc_hsdt80, hcc_hsdt81, hcc_hsdt82, hcc_hsdt83, hcc_hsdt84, hcc_hsdt85, hcc_hsdt86, hcc_hsdt87, hcc_hsdt88, hcc_hsdt89, hcc_hsdt90, hcc_hsdt91, hcc_hsdt92, hcc_hsdt93, hcc_hsdt94, hcc_hsdt95, hcc_hsdt96, hcc_hsdt97, hcc_hsdt98, hcc_hsdt99, hcc_hsdt100, hcc_hsdt101, hcc_hsdt102, hcc_hsdt103, hcc_hsdt104, hcc_hsdt105, hcc_hsdt106, hcc_hsdt107, hcc_hsdt108, hcc_hsdt109, hcc_hsdt110, hcc_hsdt111, hcc_hsdt112, hcc_hsdt113, hcc_hsdt114, hcc_hsdt115, hcc_hsdt116, hcc_hsdt117, hcc_hsdt118, hcc_hsdt119, hcc_hsdt120, hcc_hsdt121, hcc_hsdt122, hcc_hsdt123, hcc_hsdt124, hcc_hsdt125, hcc_hsdt126, hcc_hsdt127, hcc_hsdt128, hcc_hsdt129, hcc_hsdt130, hcc_hsdt131, hcc_hsdt132, hcc_hsdt133, hcc_hsdt134, hcc_hsdt135, hcc_hsdt136, hcc_hsdt137, hcc_hsdt138, hcc_hsdt139, hcc_hsdt140, hcc_hsdt141, hcc_hsdt142, hcc_hsdt143, hcc_hsdt144, hcc_hsdt145, hcc_hsdt146, hcc_hsdt147, hcc_hsdt148, hcc_hsdt149, hcc_hsdt150, hcc_hsdt151, hcc_hsdt152, hcc_hsdt153, hcc_hsdt154, hcc_hsdt155, hcc_hsdt156, hcc_hsdt157, hcc_hsdt158, hcc_hsdt159, hcc_hsdt160, hcc_hsdt161, hcc_hsdt162, hcc_hsdt163, hcc_hsdt164, hcc_hsdt165, hcc_hsdt166, hcc_hsdt167, hcc_hsdt168, hcc_hsdt169, hcc_hsdt170, hcc_hsdt171, hcc_hsdt172, hcc_hsdt173, hcc_hsdt174, hcc_hsdt175, hcc_hsdt176, hcc_hsdt177, hcc_hsdt178, hcc_hsdt179, hcc_hsdt180, hcc_hsdt181, hcc_hsdt182, hcc_hsdt183, hcc_hsdt184, hcc_hsdt185, hcc_hsdt186, hcc_hsdt187, hcc_hsdt188, hcc_hsdt189, hcc_hsdt190, hcc_hsdt191, hcc_hsdt192, hcc_hsdt193, hcc_hsdt194, hcc_hsdt195, hcc_hsdt196, hcc_hsdt197, hcc_hsdt198, hcc_hsdt199, hcc_hsdt200, hcc_hsdt201, hcc_hsdt202, hcc_hsdt203, hcc_hsdt204, hcc_hsdt205, hcc_hsdt206, hcc_hsdt207, hcc_hsdt208, hcc_hsdt209, hcc_hsdt210, hcc_hsdt211, hcc_hsdt212, hcc_hsdt213, hcc_hsdt214, hcc_hsdt215, hcc_hsdt216, hcc_hsdt217, hcc_hsdt218, hcc_hsdt219, hcc_hsdt220, hcc_hsdt221, hcc_hsdt222, hcc_hsdt223, hcc_hsdt224, hcc_hsdt225, hcc_hsdt226, hcc_hsdt227, hcc_hsdt228, hcc_hsdt229, hcc_hsdt230, hcc_hsdt231, hcc_hsdt232, hcc_hsdt233, hcc_hsdt234, hcc_hsdt235, hcc_hsdt236, hcc_hsdt237, hcc_hsdt238, hcc_hsdt239, hcc_hsdt240, hcc_hsdt241, hcc_hsdt242);
icc_hsdate=min(icc_hsdt0, icc_hsdt1, icc_hsdt2, icc_hsdt3, icc_hsdt4, icc_hsdt5, icc_hsdt6, icc_hsdt7, icc_hsdt8, icc_hsdt9, icc_hsdt10, icc_hsdt11, icc_hsdt12, icc_hsdt13, icc_hsdt14, icc_hsdt15, icc_hsdt16, icc_hsdt17, icc_hsdt18, icc_hsdt19, icc_hsdt20, icc_hsdt21, icc_hsdt22, icc_hsdt23, icc_hsdt24, icc_hsdt25, icc_hsdt26, icc_hsdt27, icc_hsdt28, icc_hsdt29, icc_hsdt30, icc_hsdt31, icc_hsdt32, icc_hsdt33, icc_hsdt34, icc_hsdt35, icc_hsdt36, icc_hsdt37, icc_hsdt38, icc_hsdt39, icc_hsdt40, icc_hsdt41, icc_hsdt42, icc_hsdt43, icc_hsdt44, icc_hsdt45, icc_hsdt46, icc_hsdt47, icc_hsdt48, icc_hsdt49, icc_hsdt50, icc_hsdt51, icc_hsdt52, icc_hsdt53, icc_hsdt54, icc_hsdt55, icc_hsdt56, icc_hsdt57, icc_hsdt58, icc_hsdt59, icc_hsdt60, icc_hsdt61, icc_hsdt62, icc_hsdt63, icc_hsdt64, icc_hsdt65, icc_hsdt66, icc_hsdt67, icc_hsdt68, icc_hsdt69, icc_hsdt70, icc_hsdt71, icc_hsdt72, icc_hsdt73, icc_hsdt74, icc_hsdt75, icc_hsdt76, icc_hsdt77, icc_hsdt78, icc_hsdt79, icc_hsdt80, icc_hsdt81, icc_hsdt82, icc_hsdt83, icc_hsdt84, icc_hsdt85, icc_hsdt86, icc_hsdt87, icc_hsdt88, icc_hsdt89, icc_hsdt90, icc_hsdt91, icc_hsdt92, icc_hsdt93, icc_hsdt94, icc_hsdt95, icc_hsdt96, icc_hsdt97, icc_hsdt98, icc_hsdt99, icc_hsdt100, icc_hsdt101, icc_hsdt102, icc_hsdt103, icc_hsdt104, icc_hsdt105, icc_hsdt106, icc_hsdt107, icc_hsdt108, icc_hsdt109, icc_hsdt110, icc_hsdt111, icc_hsdt112, icc_hsdt113, icc_hsdt114, icc_hsdt115, icc_hsdt116, icc_hsdt117, icc_hsdt118, icc_hsdt119, icc_hsdt120, icc_hsdt121, icc_hsdt122, icc_hsdt123, icc_hsdt124, icc_hsdt125, icc_hsdt126, icc_hsdt127, icc_hsdt128, icc_hsdt129, icc_hsdt130, icc_hsdt131, icc_hsdt132, icc_hsdt133, icc_hsdt134, icc_hsdt135, icc_hsdt136, icc_hsdt137, icc_hsdt138, icc_hsdt139, icc_hsdt140, icc_hsdt141, icc_hsdt142, icc_hsdt143, icc_hsdt144, icc_hsdt145, icc_hsdt146, icc_hsdt147, icc_hsdt148, icc_hsdt149, icc_hsdt150, icc_hsdt151, icc_hsdt152, icc_hsdt153, icc_hsdt154, icc_hsdt155, icc_hsdt156, icc_hsdt157, icc_hsdt158, icc_hsdt159, icc_hsdt160, icc_hsdt161, icc_hsdt162, icc_hsdt163, icc_hsdt164, icc_hsdt165, icc_hsdt166, icc_hsdt167, icc_hsdt168, icc_hsdt169, icc_hsdt170, icc_hsdt171, icc_hsdt172, icc_hsdt173, icc_hsdt174, icc_hsdt175, icc_hsdt176, icc_hsdt177, icc_hsdt178, icc_hsdt179, icc_hsdt180, icc_hsdt181, icc_hsdt182, icc_hsdt183, icc_hsdt184, icc_hsdt185, icc_hsdt186, icc_hsdt187, icc_hsdt188, icc_hsdt189, icc_hsdt190, icc_hsdt191, icc_hsdt192, icc_hsdt193, icc_hsdt194, icc_hsdt195, icc_hsdt196, icc_hsdt197, icc_hsdt198, icc_hsdt199, icc_hsdt200, icc_hsdt201, icc_hsdt202, icc_hsdt203, icc_hsdt204, icc_hsdt205, icc_hsdt206, icc_hsdt207, icc_hsdt208, icc_hsdt209, icc_hsdt210, icc_hsdt211, icc_hsdt212, icc_hsdt213, icc_hsdt214, icc_hsdt215, icc_hsdt216, icc_hsdt217, icc_hsdt218, icc_hsdt219, icc_hsdt220, icc_hsdt221, icc_hsdt222, icc_hsdt223, icc_hsdt224, icc_hsdt225, icc_hsdt226, icc_hsdt227, icc_hsdt228, icc_hsdt229, icc_hsdt230, icc_hsdt231, icc_hsdt232, icc_hsdt233, icc_hsdt234, icc_hsdt235, icc_hsdt236, icc_hsdt237, icc_hsdt238, icc_hsdt239, icc_hsdt240, icc_hsdt241, icc_hsdt242);

/*update sld using cancer registry*/
if sld=1           then do; sldnew=1; sldnewdate=slddate;end;
 else if care_HCC=1 then do; sldnew=1; sldnewdate=caredate_hcc; end;
  else                    do; sldnew=0; sldnewdate=.;            end;

/*merged cirrhosis*/
crhsnew=0; crhsnewdate=.;
if crhs_cp=1 then do; crhsnew=1; crhsnewdate=crhs_cpdate; end;
 else if crhs_dc=1 then do; crhsnew=1; crhsnewdate=crhs_dcdate; end;
 
/*Merged cancer registry and hospital data for cancer incidence*/
if care_liver=1 then do; ca_liver=1; cadate_lic=caredate_lic; end;
 else if lic_hs=1 then do; ca_liver=1; cadate_lic=lic_hsdate;    end;
  else ca_liver=0;

if care_HCC=1 then do; ca_HCC=1; cadate_hcc=caredate_hcc; end;
 else if HCC_hs=1 then do; ca_HCC=1; cadate_hcc=hcc_hsdate; end;
  else ca_HCC=0;

if care_ICC=1 then do; ca_ICC=1; cadate_icc=caredate_icc;end;
 else if ICC_hs=1 then do; ca_ICC=1; cadate_icc=icc_hsdate; end;
  else ca_ICC=0;

/*Define chronic liver disease mortality from death file*/
cld=0;
if index(s_40001_0_0,"K700")>0 or
   index(s_40001_0_0,"K701")>0 or
   index(s_40001_0_0,"K703")>0 or
   index(s_40001_0_0,"K704")>0 or
   index(s_40001_0_0,"K709")>0 or
   index(s_40001_0_0,"K740")>0 or
   index(s_40001_0_0,"K743")>0 or
   index(s_40001_0_0,"K745")>0 or
   index(s_40001_0_0,"K746")>0 
then do; cld=1; end;

 /*Coding outcomes*/
 if ca_liver=1    then enddate=cadate_lic;
 else if dthdate^=.  then enddate=dthdate;
  else if lossdate^=. then enddate=lossdate;
   else                     enddate='29feb20'd;
if enddate>'29feb20'd then do; enddate='29feb20'd; ca_liver=0;end;

if ca_HCC=1    then enddate_HCC=cadate_hcc;
 else if dthdate^=.  then enddate_HCC=dthdate;
  else if lossdate^=. then enddate_HCC=lossdate;
   else                     enddate_HCC='29feb20'd;
if enddate_HCC>'29feb20'd then do; enddate_HCC='29feb20'd; ca_HCC=0;end;

if ca_ICC=1    then enddate_ICC=cadate_icc;
 else if dthdate^=.  then enddate_ICC=dthdate;
  else if lossdate^=. then enddate_ICC=lossdate;
   else                     enddate_ICC='29feb20'd;
if enddate_ICC>'29feb20'd then do; enddate_ICC='29feb20'd; ca_ICC=0;end;

if dthdate^=.       then deathdate=dthdate; 
 else if lossdate^=. then deathdate=lossdate;
  else                     deathdate='29feb20'd;
if dthdate^=. then death=1; 
 else               death=0;

if deathdate>'29feb20'd then do; deathdate='29feb20'd; death=0;end;

/*ca_liver ca_HCC ca_ICC cadate_li dgxyear dthdate lossdate ca_esoph cadate_es ca_stomach cadate_st ca_pacreas cadate_pa ca_colon cadate_co ca_rectum cadate_re ca_crc cadate_cr
     cvddeath cancerdeath alldeath nafld cirrhosis sld clddeath*/

if nafld=1 then enddate_nf=naflddate;
 else if dthdate^=.  then enddate_nf=dthdate;
  else if lossdate^=. then enddate_nf=lossdate;
   else                     enddate_nf='29feb20'd;
if enddate_nf>'29feb20'd then do; enddate_nf='29feb20'd; nafld=0;end;

if cirrhosis=1 then enddate_ch=cirrhosisdate;
 else if dthdate^=.  then enddate_ch=dthdate;
  else if lossdate^=. then enddate_ch=lossdate;
   else                     enddate_ch='29feb20'd;
if enddate_ch>'29feb20'd then do; enddate_ch='29feb20'd; cirrhosis=0;end;

if crhsnew=1 then enddate_chn=crhsnewdate;
 else if dthdate^=.  then enddate_chn=dthdate;
  else if lossdate^=. then enddate_chn=lossdate;
   else                     enddate_chn='29feb20'd;
if enddate_chn>'29feb20'd then do; enddate_chn='29feb20'd; crhsnew=0;end;

if crhs_cp=1 then enddate_cp=crhs_cpdate;
 else if dthdate^=.  then enddate_cp=dthdate;
  else if lossdate^=. then enddate_cp=lossdate;
   else                     enddate_cp='29feb20'd;
if enddate_cp>'29feb20'd then do; enddate_cp='29feb20'd; crhs_cp=0;end;

if crhs_dc=1 then enddate_dc=crhs_dcdate;
 else if dthdate^=.  then enddate_dc=dthdate;
  else if lossdate^=. then enddate_dc=lossdate;
   else                     enddate_dc='29feb20'd;
if enddate_dc>'29feb20'd then do; enddate_dc='29feb20'd; crhs_dc=0;end;

if sld=1 then enddate_sld=slddate;
 else if dthdate^=.  then enddate_sld=dthdate;
  else if lossdate^=. then enddate_sld=lossdate;
   else                     enddate_sld='29feb20'd;
if enddate_sld>'29feb20'd then do; enddate_sld='29feb20'd; sld=0;end;

if sldnew=1 then enddate_sldnew=sldnewdate;
 else if dthdate^=.  then enddate_sldnew=dthdate;
  else if lossdate^=. then enddate_sldnew=lossdate;
   else                     enddate_sldnew='29feb20'd;
if enddate_sldnew>'29feb20'd then do; enddate_sldnew='29feb20'd; sldnew=0;end;

if cld=1 then enddate_cld=dthdate;
 else if dthdate^=.  then enddate_cld=dthdate;
  else if lossdate^=. then enddate_cld=lossdate;
   else                     enddate_cld='29feb20'd;
if enddate_cld>'29feb20'd then do; enddate_cld='29feb20'd; cld=0;end;
run;

/*************************************************************************************************************************/
/*                                              Part I : Variables readin-covariates                                     */
/*************************************************************************************************************************/

proc format;
 value agegrp 0="39-<45"
              1="45-<50"
			  2="50-<55"
			  3="55-<60"
			  4="60-<65"
			  5="65+";
 value sex 0="female"
           1="male";
 value ethnic 0="whites"
              1="non-whites";
 value mvpagrp 0="<500 METs"
               1="500-<1000 METs"
			   2="1000+ METs"
			   3="Missing";
 value bmigrp 1="<18.5"
              2="18.5-24.9"
			  3="25-29.9"
			  4="30+";
 value smk 0="never smoker"
           1="former smoker"
		   2="current smoker";
 value drk 0="never drinker"
           1="former drinker"
		   2="current drinker";
 value drkdosegrp 0="never or occasionally"
                  1="1-2/week"
				  2="3-4/week"
				  3="Almost daily";
 value yesno 0="no"
             1="yes";
 value chart 0="included"
             1="un-typical diet"
             2="baseline cancers"
			 3="baseline cirrhosis/hepatitis"
             4="Extreme energy intake"
			 5="liver cancer or died before diet questionnaire returned"
			 6="missing values in covariates";
run;
data covariate;  set datloc.covariates; run;
data medication; set datloc.medication; run;
data alcoholbsl; set datloc.alcoholbaseline;
array alctype {*} n_4407_0_0 n_4418_0_0 n_4429_0_0 n_4440_0_0 n_4451_0_0 n_4462_0_0
                  n_1568_0_0 n_1578_0_0 n_1588_0_0 n_1598_0_0 n_1608_0_0 n_5364_0_0;
do i=1 to dim(alctype);
if alctype{i} in (-1,-3) then alctype{i}=.;
end;
alcoholgrm=sum(10*sum(n_4407_0_0*2,n_4418_0_0*2,n_4429_0_0*2,n_4440_0_0*1,n_4451_0_0*1,n_4462_0_0*2)/30,
               10*sum(n_1568_0_0*2,n_1578_0_0*2,n_1588_0_0*2,n_1598_0_0*1,n_1608_0_0*1,n_5364_0_0*2)/7);
run;
%getlevels(alcoholbsl, n_eid);

proc sort data=covariate;  by n_eid; run;
proc sort data=medication; by n_eid; run;
proc sort data=alcoholbsl; by n_eid; run;
data covariates; merge covariate (in=a) medication alcoholbsl; by n_eid; if a; run;

data covariates; set covariates;
/*--------------------------------------------------------------------*/
/*Age*/
age=n_21022_0_0;
if age<45 then agegrp=0;
 else if age<50 then agegrp=1;
  else if age<55 then agegrp=2;
   else if age<60 then agegrp=3;
    else if age<65 then agegrp=4;
	 else                agegrp=5;
/*--------------------------------------------------------------------*/
/*Sex*/
sex=n_31_0_0;
/*--------------------------------------------------------------------*/
/*Ethnic*/
if n_21000_0_0 in (.,-1,-3)               then ethnic=.;
 else if n_21000_0_0 in (1,1001,1002,1003) then ethnic=0;
   else                                           ethnic=1;
/*--------------------------------------------------------------------*/
/*TDI*/
tdi=n_189_0_0;
/*--------------------------------------------------------------------*/
/*Physical activity*/
mvpa=sum(n_22038_0_0,n_22039_0_0);
if mvpa=.         then mvpagrp=3;
 else if mvpa<500  then mvpagrp=0;
  else if mvpa<1000 then mvpagrp=1;
   else                   mvpagrp=2;
/*--------------------------------------------------------------------*/
/*BMI*/
bmi=n_21001_0_0; /*more missing in n_23104_0_0*/
if bmi=. then bmigrp=.;
 else if bmi<18.5 then bmigrp=1;
  else if bmi<24.9 then bmigrp=2;
   else if bmi<29.9 then bmigrp=3;
    else                  bmigrp=4;
/*--------------------------------------------------------------------*/
/*WHR*/
whr=n_48_0_0/n_49_0_0;
wc=n_48_0_0;
/*--------------------------------------------------------------------*/
/*Medical history of cancers at baseline*/
if n_2453_0_0=0      then cancerbase=0;
 else if n_2453_0_0=1 then cancerbase=1;
  else                      cancerbase=.;
/*--------------------------------------------------------------------*/
/*HBV/HCV infection: data with HBV/HCV have no 24-hDR infromation*/
/*HBV=n_23060_0_0; HCV=n_23061_0_0;*/ 
/*--------------------------------------------------------------------*/
/*T2D */
if n_2443_0_0=0      then dm=0;
 else if n_2443_0_0=1 then dm=1;
  else                      dm=.;
/*--------------------------------------------------------------------*/
/*Smoking status*/
if n_20116_0_0=0      then smk=0;
 else if n_20116_0_0=1 then smk=1; /*Pack-years missing for former smoker: n=27660*/
  else if n_20116_0_0=2 then smk=2;
   else                       smk=.;
/*Pack-years of smoking*/
if smk=0 then smkdose=0; else smkdose=n_20161_0_0; 
/*--------------------------------------------------------------------*/
/*Aspirin use*/
asp=0; if n_6154_0_0=1 then asp=1;
/*--------------------------------------------------------------------*/
/*Alcohol*/
if n_20117_0_0 in (-3) then drk=.;
 else if n_20117_0_0=0 then drk=0;
  else if n_20117_0_0=1 then drk=1;
   else                       drk=2;
/*alcohol intake*/
drkdose=n_1558_0_0;
if drkdose in (-1,-3) then drkdose=.;

if drkdose=. then drkdosegrp=.;
 else if drkdose in (4,5,6) then drkdosegrp=0;
  else if drkdose in (3)     then drkdosegrp=1;
   else if drkdose in (2)     then drkdosegrp=2;
    else if drkdose in (1)     then drkdosegrp=3;
/*heavy drinkers*/
if drkdose=.         then alcgrm=.;
 else if alcoholgrm=. then alcgrm=0;
  else                      alcgrm=alcoholgrm;
heavydrk=0;
if sex=0 then do; if alcoholgrm>=60 then heavydrk=2; else if alcoholgrm>=30 then heavydrk=1; end;
if sex=1 then do; if alcoholgrm>=50 then heavydrk=2; else if alcoholgrm>=20 then heavydrk=1; end;

/*--------------------------------------------------------------------*/
/*Energy--will get from the diet*/

/*--------------------------------------------------------------------*/

/*Baseline liver diseases*/
array baseICD {*} n_20002_0_0 n_20002_0_1 n_20002_0_2 n_20002_0_3 n_20002_0_4 n_20002_0_5 n_20002_0_6 n_20002_0_7 n_20002_0_8 
                  n_20002_0_9 n_20002_0_10 n_20002_0_11 n_20002_0_12 n_20002_0_13 n_20002_0_14 n_20002_0_15 n_20002_0_16 
                  n_20002_0_17 n_20002_0_18 n_20002_0_19 n_20002_0_20 n_20002_0_21 n_20002_0_22 n_20002_0_23 n_20002_0_24 
                  n_20002_0_25 n_20002_0_26 n_20002_0_27 n_20002_0_28
                  /*n_20002_1_0 n_20002_1_1 n_20002_1_2 n_20002_1_3 n_20002_1_4 n_20002_1_5 n_20002_1_6 n_20002_1_7 n_20002_1_8 
                  n_20002_1_9 n_20002_1_10 n_20002_1_11 n_20002_1_12 n_20002_1_13 n_20002_1_14 n_20002_1_15*/;

balcliverd=0; bhepatitis=0; bcirrhosis=0; bcholang=0;
do i=1 to dim(baseICD);
if baseICD{i} in ("1672") then balcliverd=1;
if baseICD{i} in ("1176","1177","1178","1642","1643","1644","1645","1646") then bhepatitis=1;
if baseICD{i} in ("1179","1569") then bcirrhosis=1;
/*if baseICD{i} in ("1180","1181","1182","1183","1184","1538") then bcholang=1;*/
end;

array basedisease {*} bhepatitis balcliverd bcirrhosis /*bcholang*/;
bliverhstry=0;
do i=1 to dim(basedisease);
if basedisease{i}=1 then bliverhstry=1;
end;

array missvar {*} agegrp sex ethnic bmigrp cancerbase smk asp drk dm tdi;
ind_cova=0;
do i=1 to dim(missvar);
if missvar{i}=. then ind_cova=1;
end;
format agegrp agegrp. sex sex. ethnic ethnic. mvpagrp mvpagrp. bmigrp bmigrp. cancerbase yesno. /*HBV yesno. HCV yesno.*/ smk smk. drk drk. drkdosegrp drkdosegrp. asp yesno. dm yesno.;
keep n_eid age agegrp sex ethnic bmigrp bmi cancerbase smk asp drk whr dm mvpa mvpagrp smkdose drkdose drkdosegrp alcgrm heavydrk wc ind_cova bliverhstry bhepatitis tdi;
run;
proc means n nmiss min p25 median p75 max mean std data=covariates; var alcgrm wc; run;
proc freq data=covariates; table bhepatitis heavydrk; run;

/*************************************************************************************************************************/
/*                                              Part I : Variables readin-biomarker                                      */
/*************************************************************************************************************************/
data bio; set datloc.biomarker; 
hba1c=n_30750_0_0;    *mmol/mol;
alp=n_30610_0_0;      *U/L;
alt=n_30620_0_0;      *U/L;
ast=n_30650_0_0;      *U/L;
chol=n_30690_0_0;     *mmol/L;
ggt=n_30730_0_0;      *U/L;
crp=n_30710_0_0;      *mg/L;
glucose=n_30740_0_0;  *mmol/L;
hdl=n_30760_0_0;      *mmol/L;
ldl=n_30780_0_0;      *mmol/L;
tag=n_30870_0_0*88.57;*mmol/L to mg/dL;
tbil=n_30840_0_0;     *umol/L;
tprotein=n_30860_0_0; *g/L;
alb=n_30600_0_0*0.1;  *g/L to g/dL;
platelet=n_30080_0_0; *10^9/L;

crp_high=(crp>=3);                               *mg/L;
glucose_dsnorm=(glucose>5.6 or glucose<3.9);     *mmol/L;
hba1c_high=(hba1c>=36);                          *mmol/mol;
alp_dsnorm=(alp>129 or alp<40);                  *U/L;
alt_dsnorm=(alt>55 or alt<7);                    *U/L;
ast_dsnorm=(ast>48 or ast<8);                    *U/L;
ggt_dsnorm=(ggt>61 or ggt<8);                    *U/L;
chol_high=(chol>=6.21);                          *mmol/L;
hdl_low=(hdl<=1.03);                             *mmol/L;
ldl_high=(ldl>=4.14);                            *mmol/L;
tag_high=(tag>=200);                             *mg/dL;
tbil_dsnorm=(tbil>17 or tbil<5.1);               *umol/L;
tprotein_dsnorm=(tprotein>83 or tprotein<60);    *g/L;
alb_dsnorm=(alb>5.4 or alb<3.4);                 *g/dL;
platelet_dsnorm=(platelet>=450 or platelet<=150);*10^9/L;
/*Reference for cholesterol:
https://www.hopkinsmedicine.org/health/treatment-tests-and-therapies/lipid-panel#:~:text=Here%20are%20the%20ranges%20for,or%20above%20240%20mg%2FdL
*/
run;

data blp; set datloc.bloodpressure; run;

/*************************************************************************************************************************/
/*                                              Part II : Define study population and preparation                        */
/*************************************************************************************************************************/
%getlevels(covariates, n_eid);
%getlevels(outcome, n_eid);
%getlevels(hpylori, n_eid);
%getlevels(bio, n_eid);
%getlevels(blp, n_eid);

proc sort data=covariates; by n_eid; run; 
proc sort data=outcome;    by n_eid; run;
proc sort data=hpylori;    by n_eid; run;
proc sort data=bio;        by n_eid; run;
proc sort data=blp;        by n_eid; run;

data allinone; merge covariates outcome hpylori (in=a) bio blp; by n_eid; if a;
avg_sbp=mean(n_93_0_0, n_93_0_1);
avg_dbp=mean(n_94_0_0, n_94_0_1);

deathtimet=deathdate  -starttdate;
persontimet=enddate   -starttdate;
pt_HCC   =enddate_HCC -starttdate;
pt_ICC   =enddate_ICC -starttdate;

pt_nafld   =enddate_nf  -starttdate;
pt_crhs    =enddate_ch  -starttdate;
pt_crhsnew =enddate_chn -starttdate;
pt_sld     =enddate_sld -starttdate;
pt_sldnew  =enddate_sldnew-starttdate;
pt_death   =enddate_cld -starttdate;
pt_crhscp  =enddate_cp  -starttdate;
pt_crhsdc  =enddate_dc  -starttdate;

ind_outc=0;
array pttime {*} persontimet pt_nafld pt_crhsnew pt_crhscp pt_crhsdc pt_sldnew;
 do i=1 to dim(pttime);
  if pttime{i}<0 then ind_outc=1; 
 end; /*outcome of interests before the diet questionnaire returned*/

ind_expo=0;
if hpylori2=. then ind_expo=1;

ind_chart=0;
 if starttdate=. then ind_chart=1; 
  else if cancerbase=1 then ind_chart=2;
   else if bliverhstry=1 then ind_chart=3;
   else if ind_expo=1     then ind_chart=4;
    else if ind_outc=1     then ind_chart=5;
     else if ind_cova=1     then ind_chart=6;
	 
format ind_chart chart.;

run;

proc freq data=allinone; table ind_chart; run;

proc means n nmiss data=allinone; var starttdate; run;

proc freq data=allinone; tables ind_chart ca_liver nafld crhsnew sldnew ca_HCC ca_ICC cld/missing; run; /*follow chart*/
proc freq data=allinone; tables ca_liver nafld crhsnew sldnew ca_HCC ca_ICC cld/missing; where ind_chart=0; run; /*follow chart*/

proc print data=allinone; var ca_HCC sldnew enddate_HCC enddate_sldnew enddate_sld slddate sldnewdate starttdate persontimet pt_sldnew pt_sld;
 where .<pt_sldnew<0 and ca_HCC=1;
 format enddate_HCC mmddyy. enddate_sldnew mmddyy. enddate_sld mmddyy. sldnewdate mmddyy. slddate mmddyy. starttdate mmddyy. dthdate mmddyy.;
run;

proc freq data=allinone; table (ca_HCC ca_ICC)*ind_outc/norow nocol nopercent; run;

proc freq data=allinone; table ca_HCC ca_ICC; where ind_outc=0;run;
proc freq data=allinone; table ca_HCC ca_ICC; where ind_outc=1;run;

proc freq data=allinone; tables ca_liver care_liver nafld cirrhosis sld crhs_cp crhs_dc ca_HCC ca_ICC/missing; where ind_chart=0; run;

/*Explaination:
Some cases developped liver cancer before the first typical diet: n=42, adding thse number to liver cancer or died before diet questionnaire returned*/
proc means n nmiss min p1 p25 p50 p75 p99 max mean  data=allinone; 
 var persontimet pt_nafld pt_crhs pt_crhsnew pt_sldnew pt_crhscp pt_crhsdc pt_HCC pt_ICC pt_death deathtimet
     ca_liver ca_HCC ca_ICC 
     nafld cirrhosis crhsnew crhs_cp crhs_dc sld sldnew cld tdi; 
where ind_chart=0;
run;

data allinoneuse; set allinone; if ind_chart=0; run; /*Analysis sample: 169,516*/


libname datahpy "/udd/nhlon/ukbhpy";

data datahpy.hpyukb_11062024; set allinoneuse; run;

proc means n nmiss min p25 p50 p75 max mean  data=allinoneuse; 
var persontimet pt_nafld pt_crhs pt_crhsnew pt_sldnew pt_crhscp pt_crhsdc pt_HCC pt_ICC pt_death deathtimet
     ca_liver ca_HCC ca_ICC 
     nafld cirrhosis crhsnew crhs_cp crhs_dc sld sldnew cld
     age agegrp sex ethnic bmigrp bmi cancerbase smk asp drk heavydrk whr dm mvpa mvpagrp smkdose drkdose drkdosegrp wc bliverhstry bhepatitis tdi
     hpylori1 hpylori2;
run;

/*------------------------------------------------------------------------------------------------------------------------------*/
/*       Project name: HPYUKB                                                                                                   */
/*          Objective: To explore Helicobacter pylori and liver outcomes in the UKB                                             */
/*            Dataset: Hpylori, covariates, medication, outcomes, biomarkers                                                    */
/*         Programmer: Longgang Zhao (nhlon@channing.harvard.edu)                                                               */
/*   Program reviewer: Yun chen (yun.chen@yale.edu)                                                                             */
/*               Date: 12/08/2023 -- 01/30/2024                                                                                 */
/*           Exposure: Helicobacter pylori                                                                                      */
/*            Outcome: NAFLD, cirrhosis, severe liver disease, and liver death                                                  */
/*         Covariates: Age, sex, ethnicity, TDI, body mass index, diabetes, physical activity                                   */
/*------------------------------------------------------------------------------------------------------------------------------*/

/*Macros and notes*/
%macro histogram (datain,varlist);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne );%let var=%scan(&varlist.,&varnumber.);
title "Distribution of &var";proc sgplot data=&datain; histogram &var;run;
%let varnumber=%eval(&varnumber.+1); %end;
title;
%mend;
%macro anova(datain,varlist,group); 
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne );%let var=%scan(&varlist.,&varnumber.);
proc anova data=&datain; class &group; model &var = &group; run; 
%let varnumber=%eval(&varnumber.+1); %end;
%mend;
%macro cat(dat,cat);proc freq data=&dat; tables &cat/norow nocol nopercent;run;%mend;
%macro con(dat,con);proc means data=&dat n nmiss min p1 p25 p50 p75 p99 max; var &con; run;%mend;
/*get levels*/
%macro getlevels(dataset, var);
title "&dataset and &var";
proc sql; select count(distinct &var) as Levels,
count(*) as Nobs
from &dataset;
quit;
title;
%mend;

/*define per SD*/
%macro perSD(datain,varlist);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne ); %let var=%scan(&varlist.,&varnumber.);
/*standardization of exposure*/
proc means data=&datain; var &var; output out=stdfile; run;
proc sql; create table tmp as select &var as std from stdfile where _stat_="STD"; quit;
data &datain; if _n_=1 then set tmp (keep=std); set &datain; if std in (0,.) then &var._s=.; else &var._s=&var/std; run;
%let varnumber=%eval(&varnumber.+1);
%end;
%mend;
%macro perUnit(datain,varlist,units);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne ); %let var=%scan(&varlist.,&varnumber.);
/*standardization of exposure*/
data &datain; set &datain; &var._u=&var/&units; run;
%let varnumber=%eval(&varnumber.+1);
%end;
%mend;

/*assign median value for trend test*/
%macro median(datain,varlist,suffix=t);
%let varnumber=1; %do %while (%scan(&varlist.,&varnumber.) ne ); %let var=%scan(&varlist.,&varnumber.);
/*assign median value for each percentile*/
proc sort data=&datain; by &var.&suffix;run;
proc means data=&datain; var &var;class &var.&suffix; output out=medi median=median; run; 
data media; set medi;
if &var.&suffix=. then delete;
&var._m=median;
keep &var.&suffix &var._m;
rename &var.&suffix=group;
run;
data mlinear;set media; rename group=&var.&suffix;run;
data &datain; merge &datain mlinear; by &var.&suffix; run;
%let varnumber=%eval(&varnumber.+1);
%end;
%mend;

%macro checkvar (datain,dec,cov);
proc means n nmiss min p1 p25 median p75 p99 max mean std data=&datain maxdec=&dec; var &cov; run;
%mend;

proc format;
 value agegrp 0="39-<45"
              1="45-<50"
			  2="50-<55"
			  3="55-<60"
			  4="60-<65"
			  5="65+";
 value sex 0="female"
           1="male";
 value ethnic 0="whites"
              1="non-whites";
 value mvpagrp 0="<500 METs"
               1="500-<1000 METs"
			   2="1000+ METs"
			   3="Missing";
 value bmigrp 1="<18.5"
              2="18.5-24.9"
			  3="25-29.9"
			  4="30+";
 value smk 0="never smoker"
           1="former smoker"
		   2="current smoker";
 value drk 0="never drinker"
           1="former drinker"
		   2="current drinker";
 value drkdosegrp 0="never or occasionally"
                  1="1-2/week"
				  2="3-4/week"
				  3="Almost daily";
 value yesno 0="no"
             1="yes";
 value chart 0="included"
             1="un-typical diet"
             2="baseline cancers"
			 3="baseline cirrhosis/hepatitis"
             4="Extreme energy intake"
			 5="liver cancer or died before diet questionnaire returned"
			 6="missing values in covariates";
run;


libname datahpy "C:\Users\lz592\OneDrive - Yale University\NHANES\HpSLD\Code\UKB";
data allinone; set datahpy.hpyukb_11062024; run;

libname datloc "C:\Users\lz592\OneDrive - Yale University\NHANES\HpSLD\Code\UKB";
data pdffdata; set datloc.liverpdff; run;

%getlevels(allinone, n_eid);
%getlevels(pdffdata, n_eid);

proc sort data=allinone; by n_eid; run;
proc sort data=pdffdata; by n_eid; run;

proc means n nmiss min p25 median p75 max mean std data=allinone; var alcgrm wc bmi alcgrm heavydrk glucose avg_sbp avg_dbp tag hdl; run;

data allinonenew; merge allinone (in=a) pdffdata (in=b); 
 by n_eid; 
 if a and b;
 if n_40061_2_0 =.       then pdffnafld=.;
  else if n_40061_2_0 <5  then pdffnafld=0;
   else if n_40061_2_0 >=5 then pdffnafld=1;
run;

proc freq data=allinonenew; table pdffnafld heavydrk nafld; run;

proc means n nmiss min p25 p50 p75 max mean  data=allinonenew; 
var persontimet pt_nafld pt_crhs pt_crhsnew pt_sldnew pt_crhscp pt_crhsdc pt_HCC pt_ICC pt_death deathtimet
     ca_liver ca_HCC ca_ICC 
     nafld cirrhosis crhsnew crhs_cp crhs_dc sld sldnew cld
     age agegrp sex ethnic bmigrp bmi cancerbase smk asp drk whr dm mvpa mvpagrp smkdose drkdose drkdosegrp wc bliverhstry bhepatitis tdi;
run;

proc freq data=allinonenew; 
 table ca_liver ca_HCC ca_ICC 
       nafld cirrhosis crhsnew crhs_cp crhs_dc sld sldnew cld pdffnafld 
       hpylori1 hpylori2
       hpylori1*(ca_liver ca_HCC ca_ICC nafld cirrhosis crhsnew crhs_cp crhs_dc sld sldnew cld pdffnafld)
       hpylori2*(ca_liver ca_HCC ca_ICC nafld cirrhosis crhsnew crhs_cp crhs_dc sld sldnew cld pdffnafld);
run;

/*********************************************************************************************************/
/*Table I*/
%include "C:\Users\lz592\OneDrive - Yale University\Tools\SASMacro\Macro_describe.sas";

/*data tabledata; set allinonenew; if hpylori2 ^= . and masld in (0,1) and heavydrk=0; run;*/

proc freq data=allinonenew; table hpylori1; run;

data tabledata; set allinonenew;
if nafld=1 then newnafld=1;
 else if pdffnafld=1 then newnafld=1;
  else newnafld=0;
ind_hp=0;
if CagA=. or VacA=. or OMP=. or GroEL=. or Cata=. or UreA=. or starttdate<0 then ind_hp=1;

ind_pdff=0;
if pdffnafld =. then ind_pdff=1;
ind_all=0;

ind_hppdff=0;
if ind_pdff = 1 or ind_hp = 1 then ind_hppdff=1;

if hpylori1=0 and CagA=0 then hpcag=0;
 else if hpylori1=1 and CagA=0 then hpcag=1;
  else if hpylori1=1 and CagA=1 then hpcag=2;
   else if hpylori1=0 and CagA=1 then hpcag=.;
run;

proc means n nmiss data=tabledata; var CagA VacA OMP GroEL Cata UreA; run;

proc freq data=tabledata; table ind_all ind_hp ind_pdff ind_hppdff; run;

/*compare difference between excluded and included population*/
%describe(inputdata=tabledata,exposure=ind_all,method=0,
            varlist=age tdi bmi wc alcgrm,
            catlist=agegrp sex ethnic bmigrp smk drkdosegrp dm mvpagrp);
proc freq data=tabledata; table sex ethnic smk drkdosegrp dm mvpagrp; run;

%describe(inputdata=tabledata,exposure=ind_hp,method=0,
            varlist=age tdi bmi wc alcgrm,
            catlist=agegrp sex ethnic bmigrp smk drkdosegrp dm mvpagrp);
%describe(inputdata=tabledata,exposure=ind_hppdff,method=0,
            varlist=age tdi bmi wc alcgrm,
            catlist=agegrp sex ethnic bmigrp smk drkdosegrp dm mvpagrp);

proc freq data=tabledata; table CagA; run;

data tablecomp; set tabledata; if ind_hp=0; run;

/*definition II did not use CagA: hpylori2; so we used hpylori1*/

%describe(inputdata=tablecomp,exposure=hpylori1,method=1,
            varlist=age tdi bmi wc alcgrm,
            catlist=agegrp sex ethnic bmigrp smk drkdosegrp dm mvpagrp);

%describe(inputdata=tablecomp,exposure=CagA,method=1,
            varlist=age tdi bmi wc alcgrm,
            catlist=agegrp sex ethnic bmigrp smk drkdosegrp dm mvpagrp);

%describe(inputdata=tablecomp,exposure=hpcag,method=1,
            varlist=age tdi bmi wc alcgrm,
            catlist=agegrp sex ethnic bmigrp smk drkdosegrp dm mvpagrp);
/*total*/
proc means data=tablecomp; var age tdi bmi wc alcgrm; run;
proc freq data=tablecomp; table sex ethnic smk drkdosegrp dm mvpagrp nafld hpylori1 hpylori2; run;

proc freq data=tabledata; table sex ethnic smk drkdosegrp dm mvpagrp nafld hpylori1 hpylori2; run;

/*********************************************************************************************************/
/*Table II*/
proc freq data=tabledata; table agegrp sex ethnic pdffnafld crhsnew masld nafld newnafld; run;
proc freq data=tablecomp; table agegrp sex ethnic pdffnafld crhsnew masld nafld newnafld; run;

proc freq data=tablecomp; /*complete dataset, with CagA, use definition I: hylori1, N=4246*/
 table hpylori1*(newnafld nafld pdffnafld masld crhsnew sldnew ca_liver)/nopercent nocol norow;
run;

/*joint analysis*/
proc freq data=tablecomp; table (hpylori1 CagA hpcag)*(nafld pdffnafld newnafld)/missing norow nopercent nocol; run;

/*Table 2. Poisson to get RR*/
%macro models (datain=,exposure=,outcome=);
proc genmod data = &datain;
  class &exposure (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model &outcome =&exposure/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR' &exposure 1 -1/ exp;
run;
proc genmod data = &datain;
  class &exposure (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model &outcome =&exposure age sex ethnic tdi smk mvpagrp/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR' &exposure 1 -1/ exp;
run;
%mend;

/*nafld*/
%models(datain=tablecomp,exposure=hpylori1,outcome=nafld); /*we used, with CagA, small sample, N=4246*/
%models(datain=tablecomp,exposure=CagA,outcome=nafld);

/*Joint analysis*/
%let exposure=hpcag;
%let outcome=nafld; 
proc genmod data = tablecomp;
  class &exposure (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model &outcome =&exposure/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for hp+ cag-' &exposure 1 0 -1/ exp;
  estimate 'RR for hp+ cag+' &exposure 0 1 -1/ exp;
run;
proc genmod data = tablecomp;
  class &exposure (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model &outcome =&exposure age sex ethnic tdi smk mvpagrp/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for hp+ cag-' &exposure 1 0 -1/ exp;
  estimate 'RR for hp+ cag+' &exposure 0 1 -1/ exp;
run;

proc freq data=tablecomp; /*complete dataset, with CagA, use definition I: hylori1, N=4250*/
 table (hpylori1 CagA hpcag)*(newnafld nafld)/nopercent nocol norow;
run;

/*interaction*/
%let outcome=nafld; 
proc genmod data = tablecomp; /*P interaction=0.56*/
  class hpylori1 (ref='0') CagA (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model &outcome =hpylori1*CagA hpylori1 CagA/dist = poisson link = log type3; 
  repeated subject = n_eid/ type = unstr;
run;

proc genmod data = tablecomp; /*P interaction=0.58*/
  class hpylori1 (ref='0') CagA (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model &outcome =hpylori1*CagA hpylori1 CagA age sex ethnic tdi smk mvpagrp/dist = poisson link = log type3; 
  repeated subject = n_eid/ type = unstr;
run;


/*PDFF sample analysis-supplementary tables*/
proc freq data=tablecomp; 
 table pdffnafld;
run;

data testpdff; set tablecomp; if pdffnafld ^=.; run;
proc freq data=testpdff; table (hpylori1 CagA hpcag)*pdffnafld/nopercent nocol; run;

/*crude model*/
proc genmod data = testpdff;
  class hpylori1 (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld =hpylori1/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' hpylori1 1 -1/ exp;
run;

proc genmod data = testpdff;
  class CagA (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld =CagA/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' CagA 1 -1/ exp;
run;


proc genmod data = testpdff;
  class hpcag (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld =hpcag/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for hp+ cag-' hpcag 1 0 -1/ exp;
  estimate 'RR for hp+ cag+' hpcag 0 1 -1/ exp;
run;

/*adjusted model*/
proc genmod data = testpdff;
  class hpylori1 (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld =hpylori1 age sex ethnic tdi smk mvpagrp/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' hpylori1 1 -1/ exp;
run;

proc genmod data = testpdff;
  class CagA (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld =CagA age sex ethnic tdi smk mvpagrp/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' CagA 1 -1/ exp;
run;


proc genmod data = testpdff;
  class hpcag (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld =hpcag age sex ethnic tdi smk mvpagrp/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for hp+ cag-' hpcag 1 0 -1/ exp;
  estimate 'RR for hp+ cag+' hpcag 0 1 -1/ exp;
run;

proc freq data=tablecomp;
 table hpylori1*(CagA VacA OMP GroEL Cata UreA)/nopercent norow;
run;

/*N=4246*/
proc freq data=tablecomp;
 table (CagA VacA OMP GroEL Cata UreA)*nafld/nopercent norow nocol;
run;

%models(datain=tablecomp,exposure=CagA, outcome=nafld);
%models(datain=tablecomp,exposure=VacA, outcome=nafld);
%models(datain=tablecomp,exposure=OMP,  outcome=nafld);
%models(datain=tablecomp,exposure=GroEL,outcome=nafld);
%models(datain=tablecomp,exposure=Cata, outcome=nafld);
%models(datain=tablecomp,exposure=UreA, outcome=nafld);

/*all data, n=8484*/
proc freq data=tabledata;
 table (VacA OMP GroEL Cata UreA)*nafld/nopercent norow nocol;
run;

%models(datain=tabledata,exposure=VacA, outcome=nafld);
%models(datain=tabledata,exposure=OMP,  outcome=nafld);
%models(datain=tabledata,exposure=GroEL,outcome=nafld);
%models(datain=tabledata,exposure=Cata, outcome=nafld);
%models(datain=tabledata,exposure=UreA, outcome=nafld);


/*******************************************************************************************************************/
/*           New analysis added 2028/03/01                                                                         */
/*******************************************************************************************************************/

/*Table 4-Survival Analysis-20260225 with all new requested analyses*/
proc print data=tablecomp;
 var n_eid nafld pt_nafld enddate_nf  starttdate naflddate dthdate lossdate s_3166_0_0;
 where pt_nafld>5000;
run;

data surv; set tablecomp;
 if starttdate<0 then delete;
 pynafld=pt_nafld/365.25; 
run;

proc means data=surv n nmiss min p1 p25 p50 p75 p99 max mean std; var pynafld; run;

proc freq data=surv; table hpylori1 hpylori1*nafld; run;

/*survival curve*/
proc format;
  value hpfmt 0="H. pylori-" 1="H. pylori+";
run;

proc means data=surv min p1 p99 max; var pt_nafld; run;

ods graphics on;

/* 1) KM estimates by hpylori1 (no LIFETEST plot needed) */
proc lifetest data=surv method=km plots=none;
  time pynafld*nafld(0);
  strata hpylori1;
  ods output ProductLimitEstimates=kmplot;
run;

/* Create censor marker y only at censored time points */
data kmplot2;
  set kmplot;
  if Censor=1 then CensSurv = Survival;
  else CensSurv = .;
run;

/* Sort is recommended for step plots */
proc sort data=kmplot2;
  by Stratum pt_nafld;
run;

proc sgplot data=kmplot2;
  step    x=pynafld y=Survival / group=Stratum name="s";
  scatter x=pynafld y=CensSurv / group=Stratum markerattrs=(symbol=plus);
  yaxis min=0.98 max=1.0 label="Survival Probability";
  xaxis min=0    max=12  label="Follow-up Time (years)";
  keylegend "s" / title="H. pylori";
run;

ods graphics off;

proc freq data=surv; table death; run;
proc means data=surv; var deathtimet; run;


/**************************************************************************/
/********Cox model and different adjustments: A6, B1, B2, B4 **************/
/**************************************************************************/ 

proc freq data=surv; table bmigrp dm; run;

%macro coxmodels(datain=, exposure=, outcome=, time=pt_nafld);
ods output HazardRatios=hr1;
proc phreg data=&datain;
  class &exposure (ref='0')
        sex       (ref="female")
        ethnic    (ref="whites")
        smk       (ref="never smoker")
        mvpagrp   (ref="<500 METs")
	    drkdosegrp(ref="never or occasionally")
        / param=ref;
  model &time*&outcome(0) =
        &exposure
        age sex ethnic tdi smk mvpagrp
        / ties=efron;
  hazardratio &exposure / diff=ref cl=wald;
run;
ods output close;
data hr1fmt; set hr1;
Exposure=Description;
HRCI = cat(round(HazardRatio,0.01),' (',round(WaldLower,0.01),'-', round(WaldUpper,0.01),')');
keep Exposure HRCI;
run;
ods output HazardRatios=hr2;
proc phreg data=&datain;
  class &exposure (ref='0')
        sex       (ref="female")
        ethnic    (ref="whites")
        smk       (ref="never smoker")
        mvpagrp   (ref="<500 METs")
		drkdosegrp(ref="never or occasionally")
        / param=ref;
  model &time*&outcome(0) =
        &exposure
        age sex ethnic tdi smk mvpagrp drkdosegrp
        / ties=efron;
  hazardratio &exposure / diff=ref cl=wald;
run;
ods output close;
data hr2fmt; set hr2;
Exposure=Description;
HRCI = cat(round(HazardRatio,0.01),' (',round(WaldLower,0.01),'-', round(WaldUpper,0.01),')');
keep Exposure HRCI;
run;
ods output HazardRatios=hr3;
proc phreg data=&datain;
  class &exposure (ref='0')
        sex      (ref="female")
        ethnic   (ref="whites")
        smk      (ref="never smoker")
        mvpagrp  (ref="<500 METs")
		drkdosegrp(ref="never or occasionally")
        dm       (ref="no")
        / param=ref;
  model &time*&outcome(0) =
        &exposure
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;
  hazardratio &exposure / diff=ref cl=wald;
run;
ods output close;
data hr3fmt; set hr3;
Exposure=Description;
HRCI = cat(round(HazardRatio,0.01),' (',round(WaldLower,0.01),'-', round(WaldUpper,0.01),')');
keep Exposure HRCI;
run;

data hrfmt; set hr1fmt hr2fmt hr3fmt; run;
proc print data=hrfmt; run;
%mend;

/*nafld*/
proc freq data=surv; table (hpylori1 CagA hpcag)*nafld/nocol nopercent; run;

%coxmodels(datain=surv,exposure=hpylori1,outcome=nafld); /*N=4246*/
%coxmodels(datain=surv,exposure=CagA,outcome=nafld);
%coxmodels(datain=surv,exposure=hpcag,outcome=nafld);

/*interaction*/
proc phreg data=surv; /*P interaction=0.60*/
  class hpylori1(ref='0') CagA(ref='0') sex(ref="female") ethnic(ref="whites") smk(ref="never smoker")
        mvpagrp(ref="<500 METs") drkdosegrp(ref="never or occasionally") dm(ref="no")
        / param=ref;
  model pt_nafld*nafld(0) = hpylori1*CagA hpylori1 CagA/ ties=efron;
run;
proc phreg data=surv; /*P interaction=0.66*/
  class hpylori1(ref='0') CagA(ref='0') sex(ref="female") ethnic(ref="whites") smk(ref="never smoker")
        mvpagrp(ref="<500 METs") drkdosegrp(ref="never or occasionally") dm(ref="no")
        / param=ref;
  model pt_nafld*nafld(0) = hpylori1*CagA hpylori1 CagA age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi / ties=efron;
run;

/*adjust WC*/
proc phreg data=surv;
  class hpylori1(ref='0') sex(ref="female") ethnic(ref="whites") smk(ref="never smoker")
        mvpagrp(ref="<500 METs") drkdosegrp(ref="never or occasionally") dm(ref="no")
        / param=ref;
  model pt_nafld*nafld(0) = hpylori1 age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi wc/ ties=efron;
  hazardratio hpylori1 / diff=ref cl=wald;
run;
proc phreg data=surv;
  class CagA(ref='0') sex(ref="female") ethnic(ref="whites") smk(ref="never smoker")
        mvpagrp(ref="<500 METs") drkdosegrp(ref="never or occasionally") dm(ref="no")
        / param=ref;
  model pt_nafld*nafld(0) = CagA age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi wc/ ties=efron;
  hazardratio CagA / diff=ref cl=wald;
run;
proc phreg data=surv;
  class hpcag(ref='0') sex(ref="female") ethnic(ref="whites") smk(ref="never smoker")
        mvpagrp(ref="<500 METs") drkdosegrp(ref="never or occasionally") dm(ref="no")
        / param=ref;
  model pt_nafld*nafld(0) = hpcag age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi wc/ ties=efron;
  hazardratio hpcag / diff=ref cl=wald;
run;

/*Hp based on 5 anti-agents*/
proc freq data=surv; table hpylori1*hpylori2; run;
%coxmodels(datain=surv,exposure=hpylori2,outcome=nafld); /*N=4246*/

/*Other sensitivity tables*/
/*Other anti-agents*/
proc freq data=surv;
 table (VacA OMP GroEL Cata UreA)*nafld/nopercent norow nocol;
run;

%coxmodels(datain=surv,exposure=VacA,outcome=nafld);
%coxmodels(datain=surv,exposure=OMP,outcome=nafld);
%coxmodels(datain=surv,exposure=GroEL,outcome=nafld);
%coxmodels(datain=surv,exposure=Cata,outcome=nafld);
%coxmodels(datain=surv,exposure=UreA,outcome=nafld);

/**************************************************************************/
/*Lag analyses: B2/3 */
data lag; set surv; flowupyear=pt_nafld/365.25; run;

proc univariate data=lag; var flowupyear; histogram flowupyear; run;

proc univariate data=lag; var flowupyear; histogram flowupyear; where nafld=1; run;

data lag2; set lag; if flowupyear<2 and nafld=1 then delete; run; proc freq data=lag2; table nafld; run;
proc freq data=lag2; table (hpylori1 CagA hpcag)*nafld/norow nocol nopercent; run;

%coxmodels(datain=lag2,exposure=hpylori1,outcome=nafld);
%coxmodels(datain=lag2,exposure=CagA,outcome=nafld);
%coxmodels(datain=lag2,exposure=hpcag,outcome=nafld);

data lag3; set lag; if flowupyear<3 and nafld=1 then delete; run; proc freq data=lag2; table nafld; run;
proc freq data=lag3; table (hpylori1 CagA hpcag)*nafld/norow nocol nopercent; run;

%coxmodels(datain=lag3,exposure=hpylori1,outcome=nafld);
%coxmodels(datain=lag3,exposure=CagA,outcome=nafld);
%coxmodels(datain=lag3,exposure=hpcag,outcome=nafld);


/**************************************************************************/
/*Landmark analyses: B5 */
data lag5; set lag; if flowupyear<5 and nafld=1 then delete; run; proc freq data=lag5; table nafld; run;
proc freq data=lag5; table (hpylori1 CagA hpcag)*nafld/norow nocol nopercent; run;

%coxmodels(datain=lag5,exposure=hpylori1,outcome=nafld);
%coxmodels(datain=lag5,exposure=CagA,outcome=nafld);
%coxmodels(datain=lag5,exposure=hpcag,outcome=nafld);

/**************************************************************************/
/********Charicteristics between population: B6              **************/
/**************************************************************************/ 
*data showed above with the baseline table 1;

/**************************************************************************/
/********IPTW analysis: B6                                   **************/
/**************************************************************************/

/*******IPTW for serology sample  ************************/
/*========================================================
  0) Define inclusion flag: ind_hp=0 means included
========================================================*/
data tabledata_ipw;
  set tabledata;
  incl = (ind_hp = 0);   /* 1=in analytic sample, 0=excluded */
run;

/*========================================================
  1) Inclusion model within serology sample (N=8,448)
========================================================*/
proc logistic data=tabledata_ipw descending;
  class sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp(ref="never or occasionally")
        dm         (ref="no")
        / param=ref;
  model incl(event='1') =
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi;
  output out=ipw_pred pred=ps_incl;
run;

/* Stabilization numerator: overall inclusion probability */
proc means data=ipw_pred noprint;
  var incl;
  output out=_pbar mean=pbar;
run;

/* Stabilized weights for included participants only */
data ipw_wt;
  if _n_=1 then set _pbar;
  set ipw_pred;

  if incl=1 then sw = pbar / ps_incl;
  else sw = .;
run;

/* Optional truncation at 1st/99th percentiles */
proc univariate data=ipw_wt noprint;
  where incl=1 and sw ne .;
  var sw;
  output out=_pct pctlpts=1 99 pctlpre=p_;
run;

data ipw_wt;
  if _n_=1 then set _pct;
  set ipw_wt;

  sw_trunc = sw;
  if incl=1 then do;
    if sw_trunc < p_1  then sw_trunc = p_1;
    if sw_trunc > p_99 then sw_trunc = p_99;
  end;
run;

/* Weight diagnostics */
proc means data=ipw_wt n mean std min p1 p5 p50 p95 p99 max;
  where incl=1 and sw_trunc ne .;
  var sw sw_trunc;
run;

/*========================================================
  2) Weighted Cox model in analytic sample (incl=1)
========================================================*/
/*old for reference*/
proc phreg data=ipw_wt(where=(incl=1));
  class hpylori1   (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp(ref="never or occasionally")
        dm         (ref="no")
        / param=ref;
  model pt_nafld*nafld(0) =
        hpylori1
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;
  hazardratio hpylori1 / diff=ref cl=wald;
run;

%macro iptw_phreg(exposure=, ref=0, outHR=);
  %if %length(&outHR) %then %do;
    ods output HazardRatios=&outHR;
  %end;

  proc phreg data=ipw_wt(where=(incl=1));
    class &exposure (ref="&ref")
          sex        (ref="female")
          ethnic     (ref="whites")
          smk        (ref="never smoker")
          mvpagrp    (ref="<500 METs")
          drkdosegrp (ref="never or occasionally")
          dm         (ref="no")
          / param=ref;

    weight sw_trunc;
    id n_eid;

    model pt_nafld*nafld(0) =
          &exposure
          age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
          / ties=efron;
    hazardratio &exposure / diff=ref cl=wald;
  run;

  %if %length(&outHR) %then %do;
    ods output close;
  %end;
%mend;

/* Run 3 models */
%iptw_phreg(exposure=hpylori1, ref=0, outHR=HR_hp_w);
%iptw_phreg(exposure=CagA,     ref=0, outHR=HR_caga_w);
%iptw_phreg(exposure=hpcag,    ref=0, outHR=HR_joint_w);


/*******IPTW for PDFF sample  ************************/
/*========================================================
  0) Define inclusion flag: ind_hp=0 means included
========================================================*/
proc freq data=tabledata; table ind_hppdff*pdffnafld; run;

data tabledata_ipw;
  set tabledata;
  incl = (ind_hppdff = 0);   /* 1=in analytic sample, 0=excluded */
run;

/*========================================================
  1) Inclusion model within serology sample (N=8,448)
========================================================*/
proc logistic data=tabledata_ipw descending;
  class sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp(ref="never or occasionally")
        dm         (ref="no")
        / param=ref;
  model incl(event='1') =
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi;
  output out=ipw_pred pred=ps_incl;
run;

/* Stabilization numerator: overall inclusion probability */
proc means data=ipw_pred noprint;
  var incl;
  output out=_pbar mean=pbar;
run;

/* Stabilized weights for included participants only */
data ipw_wt;
  if _n_=1 then set _pbar;
  set ipw_pred;

  if incl=1 then sw = pbar / ps_incl;
  else sw = .;
run;

/* Optional truncation at 1st/99th percentiles */
proc univariate data=ipw_wt noprint;
  where incl=1 and sw ne .;
  var sw;
  output out=_pct pctlpts=1 99 pctlpre=p_;
run;

data ipw_wt;
  if _n_=1 then set _pct;
  set ipw_wt;

  sw_trunc = sw;
  if incl=1 then do;
    if sw_trunc < p_1  then sw_trunc = p_1;
    if sw_trunc > p_99 then sw_trunc = p_99;
  end;
run;

/* Weight diagnostics */
proc means data=ipw_wt n mean std min p1 p5 p50 p95 p99 max;
  where incl=1 and sw_trunc ne .;
  var sw sw_trunc;
run;

/*========================================================
  2) Weighted Cox model in analytic sample (incl=1)
========================================================*/
/*old for reference*/
proc freq data=testpdff; table (hpylori1 CagA hpcag)*pdffnafld/norow nocol nopercent; run;
proc genmod data = testpdff;
  class hpylori1 (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld=hpylori1 age sex ethnic tdi smk mvpagrp drkdosegrp/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' hpylori1 1 -1/ exp;
run;

proc genmod data = testpdff;
  class hpylori1 (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld=hpylori1 age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' hpylori1 1 -1/ exp;
run;

proc genmod data = testpdff;
  class CagA (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld=CagA age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
  estimate 'RR for + vs-' CagA 1 -1/ exp;
run;

proc genmod data = testpdff;
  class hpcag (ref='0') agegrp (ref="45-<50") sex (ref="female") ethnic (ref="whites") smk (ref="never smoker") mvpagrp (ref="<500 METs") n_eid; 
  model pdffnafld=hpcag age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi/dist = poisson link = log; 
  repeated subject = n_eid/ type = unstr;
    estimate 'RR for hp+ cag-' hpcag 1 0 -1/ exp;
    estimate 'RR for hp+ cag+' hpcag 0 1 -1/ exp;
run;

%macro iptw_poisson_pdff(exposure=, ref=0, outcome=, outRR=, corr=ind);
  %if %length(&outRR) %then %do;
    ods output Estimates=&outRR;
  %end;

  proc genmod data=ipw_wt(where=(incl=1)) descending;
    class &exposure (ref="&ref")
          sex        (ref="female")
          ethnic     (ref="whites")
          smk        (ref="never smoker")
          mvpagrp    (ref="<500 METs")
          drkdosegrp (ref="never or occasionally")
          dm         (ref="no")
		  n_eid
          / param=ref;
    model &outcome =
          &exposure
          age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
          / dist=poisson link=log;

    /* IPTW weight */
    weight sw_trunc;
    /* Use GEE-style robust SEs clustered by participant */
    repeated subject=n_eid / type=&corr;

    /* RR for exposure 1 vs ref */
    estimate "RR for &exposure vs ref" &exposure 1 / exp;
  run;

  %if %length(&outRR) %then %do;
    ods output close;
  %end;

%mend;

/* Example calls */
%iptw_poisson_pdff(exposure=hpylori1, ref=0, outcome=pdffnafld, outRR=RR_hp_pdff_w);
%iptw_poisson_pdff(exposure=CagA,     ref=0, outcome=pdffnafld, outRR=RR_caga_pdff_w);

proc genmod data=ipw_wt(where=(incl=1)) descending;
    class hpcag      (ref="0")
          sex        (ref="female")
          ethnic     (ref="whites")
          smk        (ref="never smoker")
          mvpagrp    (ref="<500 METs")
          drkdosegrp (ref="never or occasionally")
          dm         (ref="no")
		  n_eid
          / param=ref;
    model pdffnafld =
          hpcag
          age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
          / dist=poisson link=log;
    weight sw_trunc;
    repeated subject=n_eid / type=ind;
    estimate 'RR for hp+ cag-' hpcag 1 0 -1/ exp;
    estimate 'RR for hp+ cag+' hpcag 0 1 -1/ exp;
  run;

/**************************************************************************/
/********Negative control analysis: B9                       **************/
/**************************************************************************/ 

/*******Negative outcome control using osteoarthritis ************************/
data negcontroloutcome; set tabledata; 
/*Define osteo from hospital data*/
array hspicd10 {*} s_41270_0_0 s_41270_0_1 s_41270_0_2 s_41270_0_3 s_41270_0_4 s_41270_0_5 s_41270_0_6 s_41270_0_7 s_41270_0_8 s_41270_0_9 s_41270_0_10 
                  s_41270_0_11 s_41270_0_12 s_41270_0_13 s_41270_0_14 s_41270_0_15 s_41270_0_16 s_41270_0_17 s_41270_0_18 s_41270_0_19 s_41270_0_20 
                  s_41270_0_21 s_41270_0_22 s_41270_0_23 s_41270_0_24 s_41270_0_25 s_41270_0_26 s_41270_0_27 s_41270_0_28 s_41270_0_29 s_41270_0_30 
                  s_41270_0_31 s_41270_0_32 s_41270_0_33 s_41270_0_34 s_41270_0_35 s_41270_0_36 s_41270_0_37 s_41270_0_38 s_41270_0_39 s_41270_0_40 
                  s_41270_0_41 s_41270_0_42 s_41270_0_43 s_41270_0_44 s_41270_0_45 s_41270_0_46 s_41270_0_47 s_41270_0_48 s_41270_0_49 s_41270_0_50 
                  s_41270_0_51 s_41270_0_52 s_41270_0_53 s_41270_0_54 s_41270_0_55 s_41270_0_56 s_41270_0_57 s_41270_0_58 s_41270_0_59 s_41270_0_60 
                  s_41270_0_61 s_41270_0_62 s_41270_0_63 s_41270_0_64 s_41270_0_65 s_41270_0_66 s_41270_0_67 s_41270_0_68 s_41270_0_69 s_41270_0_70 
                  s_41270_0_71 s_41270_0_72 s_41270_0_73 s_41270_0_74 s_41270_0_75 s_41270_0_76 s_41270_0_77 s_41270_0_78 s_41270_0_79 s_41270_0_80 
                  s_41270_0_81 s_41270_0_82 s_41270_0_83 s_41270_0_84 s_41270_0_85 s_41270_0_86 s_41270_0_87 s_41270_0_88 s_41270_0_89 s_41270_0_90 
                  s_41270_0_91 s_41270_0_92 s_41270_0_93 s_41270_0_94 s_41270_0_95 s_41270_0_96 s_41270_0_97 s_41270_0_98 s_41270_0_99 s_41270_0_100 
                  s_41270_0_101 s_41270_0_102 s_41270_0_103 s_41270_0_104 s_41270_0_105 s_41270_0_106 s_41270_0_107 s_41270_0_108 s_41270_0_109 s_41270_0_110 
                  s_41270_0_111 s_41270_0_112 s_41270_0_113 s_41270_0_114 s_41270_0_115 s_41270_0_116 s_41270_0_117 s_41270_0_118 s_41270_0_119 s_41270_0_120 
                  s_41270_0_121 s_41270_0_122 s_41270_0_123 s_41270_0_124 s_41270_0_125 s_41270_0_126 s_41270_0_127 s_41270_0_128 s_41270_0_129 s_41270_0_130 
                  s_41270_0_131 s_41270_0_132 s_41270_0_133 s_41270_0_134 s_41270_0_135 s_41270_0_136 s_41270_0_137 s_41270_0_138 s_41270_0_139 s_41270_0_140 
                  s_41270_0_141 s_41270_0_142 s_41270_0_143 s_41270_0_144 s_41270_0_145 s_41270_0_146 s_41270_0_147 s_41270_0_148 s_41270_0_149 s_41270_0_150 
                  s_41270_0_151 s_41270_0_152 s_41270_0_153 s_41270_0_154 s_41270_0_155 s_41270_0_156 s_41270_0_157 s_41270_0_158 s_41270_0_159 s_41270_0_160 
                  s_41270_0_161 s_41270_0_162 s_41270_0_163 s_41270_0_164 s_41270_0_165 s_41270_0_166 s_41270_0_167 s_41270_0_168 s_41270_0_169 s_41270_0_170 
                  s_41270_0_171 s_41270_0_172 s_41270_0_173 s_41270_0_174 s_41270_0_175 s_41270_0_176 s_41270_0_177 s_41270_0_178 s_41270_0_179 s_41270_0_180 
                  s_41270_0_181 s_41270_0_182 s_41270_0_183 s_41270_0_184 s_41270_0_185 s_41270_0_186 s_41270_0_187 s_41270_0_188 s_41270_0_189 s_41270_0_190 
                  s_41270_0_191 s_41270_0_192 s_41270_0_193 s_41270_0_194 s_41270_0_195 s_41270_0_196 s_41270_0_197 s_41270_0_198 s_41270_0_199 s_41270_0_200 
                  s_41270_0_201 s_41270_0_202 s_41270_0_203 s_41270_0_204 s_41270_0_205 s_41270_0_206 s_41270_0_207 s_41270_0_208 s_41270_0_209 s_41270_0_210 
                  s_41270_0_211 s_41270_0_212 s_41270_0_213 s_41270_0_214 s_41270_0_215 s_41270_0_216 s_41270_0_217 s_41270_0_218 s_41270_0_219 s_41270_0_220 
                  s_41270_0_221 s_41270_0_222 s_41270_0_223 s_41270_0_224 s_41270_0_225 s_41270_0_226 s_41270_0_227 s_41270_0_228 s_41270_0_229 s_41270_0_230 
                  s_41270_0_231 s_41270_0_232 s_41270_0_233 s_41270_0_234 s_41270_0_235 s_41270_0_236 s_41270_0_237 s_41270_0_238 s_41270_0_239 s_41270_0_240 
                  s_41270_0_241 s_41270_0_242;
array hspidt {*} s_41280_0_0 s_41280_0_1 s_41280_0_2 s_41280_0_3 s_41280_0_4 s_41280_0_5 s_41280_0_6 s_41280_0_7 s_41280_0_8 s_41280_0_9 s_41280_0_10 
                  s_41280_0_11 s_41280_0_12 s_41280_0_13 s_41280_0_14 s_41280_0_15 s_41280_0_16 s_41280_0_17 s_41280_0_18 s_41280_0_19 s_41280_0_20 
                  s_41280_0_21 s_41280_0_22 s_41280_0_23 s_41280_0_24 s_41280_0_25 s_41280_0_26 s_41280_0_27 s_41280_0_28 s_41280_0_29 s_41280_0_30 
                  s_41280_0_31 s_41280_0_32 s_41280_0_33 s_41280_0_34 s_41280_0_35 s_41280_0_36 s_41280_0_37 s_41280_0_38 s_41280_0_39 s_41280_0_40 
                  s_41280_0_41 s_41280_0_42 s_41280_0_43 s_41280_0_44 s_41280_0_45 s_41280_0_46 s_41280_0_47 s_41280_0_48 s_41280_0_49 s_41280_0_50 
                  s_41280_0_51 s_41280_0_52 s_41280_0_53 s_41280_0_54 s_41280_0_55 s_41280_0_56 s_41280_0_57 s_41280_0_58 s_41280_0_59 s_41280_0_60 
                  s_41280_0_61 s_41280_0_62 s_41280_0_63 s_41280_0_64 s_41280_0_65 s_41280_0_66 s_41280_0_67 s_41280_0_68 s_41280_0_69 s_41280_0_70 
                  s_41280_0_71 s_41280_0_72 s_41280_0_73 s_41280_0_74 s_41280_0_75 s_41280_0_76 s_41280_0_77 s_41280_0_78 s_41280_0_79 s_41280_0_80 
                  s_41280_0_81 s_41280_0_82 s_41280_0_83 s_41280_0_84 s_41280_0_85 s_41280_0_86 s_41280_0_87 s_41280_0_88 s_41280_0_89 s_41280_0_90 
                  s_41280_0_91 s_41280_0_92 s_41280_0_93 s_41280_0_94 s_41280_0_95 s_41280_0_96 s_41280_0_97 s_41280_0_98 s_41280_0_99 s_41280_0_100 
                  s_41280_0_101 s_41280_0_102 s_41280_0_103 s_41280_0_104 s_41280_0_105 s_41280_0_106 s_41280_0_107 s_41280_0_108 s_41280_0_109 s_41280_0_110 
                  s_41280_0_111 s_41280_0_112 s_41280_0_113 s_41280_0_114 s_41280_0_115 s_41280_0_116 s_41280_0_117 s_41280_0_118 s_41280_0_119 s_41280_0_120 
                  s_41280_0_121 s_41280_0_122 s_41280_0_123 s_41280_0_124 s_41280_0_125 s_41280_0_126 s_41280_0_127 s_41280_0_128 s_41280_0_129 s_41280_0_130 
                  s_41280_0_131 s_41280_0_132 s_41280_0_133 s_41280_0_134 s_41280_0_135 s_41280_0_136 s_41280_0_137 s_41280_0_138 s_41280_0_139 s_41280_0_140 
                  s_41280_0_141 s_41280_0_142 s_41280_0_143 s_41280_0_144 s_41280_0_145 s_41280_0_146 s_41280_0_147 s_41280_0_148 s_41280_0_149 s_41280_0_150 
                  s_41280_0_151 s_41280_0_152 s_41280_0_153 s_41280_0_154 s_41280_0_155 s_41280_0_156 s_41280_0_157 s_41280_0_158 s_41280_0_159 s_41280_0_160 
                  s_41280_0_161 s_41280_0_162 s_41280_0_163 s_41280_0_164 s_41280_0_165 s_41280_0_166 s_41280_0_167 s_41280_0_168 s_41280_0_169 s_41280_0_170 
                  s_41280_0_171 s_41280_0_172 s_41280_0_173 s_41280_0_174 s_41280_0_175 s_41280_0_176 s_41280_0_177 s_41280_0_178 s_41280_0_179 s_41280_0_180 
                  s_41280_0_181 s_41280_0_182 s_41280_0_183 s_41280_0_184 s_41280_0_185 s_41280_0_186 s_41280_0_187 s_41280_0_188 s_41280_0_189 s_41280_0_190 
                  s_41280_0_191 s_41280_0_192 s_41280_0_193 s_41280_0_194 s_41280_0_195 s_41280_0_196 s_41280_0_197 s_41280_0_198 s_41280_0_199 s_41280_0_200 
                  s_41280_0_201 s_41280_0_202 s_41280_0_203 s_41280_0_204 s_41280_0_205 s_41280_0_206 s_41280_0_207 s_41280_0_208 s_41280_0_209 s_41280_0_210 
                  s_41280_0_211 s_41280_0_212 s_41280_0_213 s_41280_0_214 s_41280_0_215 s_41280_0_216 s_41280_0_217 s_41280_0_218 s_41280_0_219 s_41280_0_220 
                  s_41280_0_221 s_41280_0_222 s_41280_0_223 s_41280_0_224 s_41280_0_225 s_41280_0_226 s_41280_0_227 s_41280_0_228 s_41280_0_229 s_41280_0_230 
                  s_41280_0_231 s_41280_0_232 s_41280_0_233 s_41280_0_234 s_41280_0_235 s_41280_0_236 s_41280_0_237 s_41280_0_238 s_41280_0_239 s_41280_0_240 
                  s_41280_0_241 s_41280_0_242;
array osteodt    {*} osteodt0-osteodt242;

osteo=0; 
do i=1 to dim(hspicd10);
 if substr(hspicd10{i},1,3) in ("M15","M16","M17","M18","M19") then do; osteo = 1; osteodt{i} = hspidt{i}; end; else osteodt{i} = .;end;
osteodate=min(osteodt0, osteodt1, osteodt2, osteodt3, osteodt4, osteodt5, osteodt6, osteodt7, osteodt8, osteodt9, osteodt10, osteodt11, osteodt12, osteodt13, osteodt14, osteodt15, osteodt16, osteodt17, osteodt18, osteodt19, osteodt20, osteodt21, osteodt22, osteodt23, osteodt24, osteodt25, osteodt26, osteodt27, osteodt28, osteodt29, osteodt30, osteodt31, osteodt32, osteodt33, osteodt34, osteodt35, osteodt36, osteodt37, osteodt38, osteodt39, osteodt40, osteodt41, osteodt42, osteodt43, osteodt44, osteodt45, osteodt46, osteodt47, osteodt48, osteodt49, osteodt50, osteodt51, osteodt52, osteodt53, osteodt54, osteodt55, osteodt56, osteodt57, osteodt58, osteodt59, osteodt60, osteodt61, osteodt62, osteodt63, osteodt64, osteodt65, osteodt66, osteodt67, osteodt68, osteodt69, osteodt70, osteodt71, osteodt72, osteodt73, osteodt74, osteodt75, osteodt76, osteodt77, osteodt78, osteodt79, osteodt80, osteodt81, osteodt82, osteodt83, osteodt84, osteodt85, osteodt86, osteodt87, osteodt88, osteodt89, osteodt90, osteodt91, osteodt92, osteodt93, osteodt94, osteodt95, osteodt96, osteodt97, osteodt98, osteodt99, osteodt100, osteodt101, osteodt102, osteodt103, osteodt104, osteodt105, osteodt106, osteodt107, osteodt108, osteodt109, osteodt110, osteodt111, osteodt112, osteodt113, osteodt114, osteodt115, osteodt116, osteodt117, osteodt118, osteodt119, osteodt120, osteodt121, osteodt122, osteodt123, osteodt124, osteodt125, osteodt126, osteodt127, osteodt128, osteodt129, osteodt130, osteodt131, osteodt132, osteodt133, osteodt134, osteodt135, osteodt136, osteodt137, osteodt138, osteodt139, osteodt140, osteodt141, osteodt142, osteodt143, osteodt144, osteodt145, osteodt146, osteodt147, osteodt148, osteodt149, osteodt150, osteodt151, osteodt152, osteodt153, osteodt154, osteodt155, osteodt156, osteodt157, osteodt158, osteodt159, osteodt160, osteodt161, osteodt162, osteodt163, osteodt164, osteodt165, osteodt166, osteodt167, osteodt168, osteodt169, osteodt170, osteodt171, osteodt172, osteodt173, osteodt174, osteodt175, osteodt176, osteodt177, osteodt178, osteodt179, osteodt180, osteodt181, osteodt182, osteodt183, osteodt184, osteodt185, osteodt186, osteodt187, osteodt188, osteodt189, osteodt190, osteodt191, osteodt192, osteodt193, osteodt194, osteodt195, osteodt196, osteodt197, osteodt198, osteodt199, osteodt200, osteodt201, osteodt202, osteodt203, osteodt204, osteodt205, osteodt206, osteodt207, osteodt208, osteodt209, osteodt210, osteodt211, osteodt212, osteodt213, osteodt214, osteodt215, osteodt216, osteodt217, osteodt218, osteodt219, osteodt220, osteodt221, osteodt222, osteodt223, osteodt224, osteodt225, osteodt226, osteodt227, osteodt228, osteodt229, osteodt230, osteodt231, osteodt232, osteodt233, osteodt234, osteodt235, osteodt236, osteodt237, osteodt238, osteodt239, osteodt240, osteodt241, osteodt242);

if osteo=1 then enddate_os=osteodate;
 else if dthdate^=.  then enddate_os=dthdate;
  else if lossdate^=. then enddate_os=lossdate;
   else                     enddate_os='29feb20'd;
if enddate_os>'29feb20'd then do; enddate_os='29feb20'd; nafld=0;end;

pt_osteo   =enddate_os  -starttdate;
run;

proc freq data=negcontroloutcome; table osteo; run;

proc means data=negcontroloutcome; var pt_osteo; run;

data negcontroloutcomeuse; set negcontroloutcome; if pt_osteo>0; run;

proc freq data=negcontroloutcomeuse; table osteo*nafld/chisq norow nopercent; run;

proc phreg data=negcontroloutcomeuse;
  class hpylori1 (ref='0')
        sex      (ref="female")
        ethnic   (ref="whites")
        smk      (ref="never smoker")
        mvpagrp  (ref="<500 METs")
		drkdosegrp(ref="never or occasionally")
        dm       (ref="no")
        / param=ref;
  model pt_osteo*osteo(0) =
        hpylori1
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;
  hazardratio hpylori1 / diff=ref cl=wald;
run;

proc phreg data=negcontroloutcomeuse;
  class CagA (ref='0')
        sex      (ref="female")
        ethnic   (ref="whites")
        smk      (ref="never smoker")
        mvpagrp  (ref="<500 METs")
		drkdosegrp(ref="never or occasionally")
        dm       (ref="no")
        / param=ref;
  model pt_osteo*osteo(0) =
        CagA
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;
  hazardratio CagA / diff=ref cl=wald;
run;

proc phreg data=negcontroloutcomeuse;
  class hpcag (ref='0')
        sex      (ref="female")
        ethnic   (ref="whites")
        smk      (ref="never smoker")
        mvpagrp  (ref="<500 METs")
		drkdosegrp(ref="never or occasionally")
        dm       (ref="no")
        / param=ref;
  model pt_osteo*osteo(0) =
        hpcag
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;
  hazardratio hpcag / diff=ref cl=wald;
run;


/*******Negative exposure control using CMV ************************/
libname antigen "C:\Users\lz592\OneDrive - Yale University\NHANES\HpSLD\Code\UKB\CEBP R1";
data antigens; set antigen.antigens; 
  cmv_pos = n_23054_0_0;
  ebv_pos  = n_23053_0_0;
  hsv1_pos = n_23050_0_0;
  array posantigen {3} ;
  do i=1 to 3;
  if posantigen{i}=. then posantigen{i}=0; 
  end;
run;

proc sort data=antigens;  by n_eid; run;
proc sort data=tabledata; by n_eid; run;

data negcontrolexposure; merge tabledata(in=a) antigens(in=b); by n_eid; if a and b; run;

proc freq data=negcontrolexposure; table hpylori1 cmv_pos ebv_pos hsv1_pos; run;

proc freq data=negcontrolexposure; table (cmv_pos ebv_pos hsv1_pos)*(hpylori1 nafld)/norow nopercent; run;

%coxmodels(datain=negcontrolexposure,exposure=cmv_pos,outcome=nafld);
%coxmodels(datain=negcontrolexposure,exposure=ebv_pos,outcome=nafld);
%coxmodels(datain=negcontrolexposure,exposure=hsv1_pos,outcome=nafld);
/*
3 cmv_pos 1 vs 0 1.59 (0.92, 2.75) 
3 ebv_pos 1 vs 0 0.9 (0.28, 2.9) 
3 hsv1_pos 1 vs 0 1.44 (0.8, 2.61) 
*/

/**************************************************************************/
/********Competing risk analysis: B10                        **************/
/**************************************************************************/ 

proc freq data=tablecomp; table nafld*death; run;
%coxmodels(datain=tablecomp,exposure=hpylori1,outcome=death,time=deathtimet);

/*--------------------------------------------------------------
  0) Create a 3-level event variable for PHREG eventcode=
     status = 0: censored
     status = 1: MASLD event (event of interest)
     status = 2: death (competing event)
--------------------------------------------------------------*/
data cr_dat;
  set tablecomp;  
  length status 8;
  if nafld=1 then status=1;
  else if death=1 then status=2;
  else status=0;
run;
proc freq data=cr_dat; table status; run;

/*--------------------------------------------------------------
  1) Cause-specific Cox model (death treated as censoring)
     (This estimates cause-specific HR for MASLD)
--------------------------------------------------------------*/
proc phreg data=cr_dat;
  class hpylori1   (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp (ref="never or occasionally")
        dm         (ref="no")
        / param=ref;
  model pt_nafld*status(0,2) =
        hpylori1
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;

  hazardratio hpylori1 / diff=ref cl=wald;
run;

proc phreg data=cr_dat;
  class CagA       (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp (ref="never or occasionally")
        dm         (ref="no")
        / param=ref;
  model pt_nafld*status(0,2) =
        CagA
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;

  hazardratio CagA / diff=ref cl=wald;
run;

proc phreg data=cr_dat;
  class hpcag      (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp (ref="never or occasionally")
        dm         (ref="no")
        / param=ref;
  model pt_nafld*status(0,2) =
        hpcag
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=efron;

  hazardratio hpcag / diff=ref cl=wald;
run;

/*--------------------------------------------------------------
  2) Fine–Gray subdistribution hazards model
     (death is a competing event)
     eventcode=1 specifies MASLD as the event of interest.
--------------------------------------------------------------*/
proc phreg data=cr_dat;
  class hpylori1   (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp (ref="never or occasionally")
        dm         (ref="no")
        / param=ref;

  model pt_nafld*status(0) =
        hpylori1
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / eventcode=1 ties=efron;

  hazardratio hpylori1 / diff=ref cl=wald;
run;

proc phreg data=cr_dat;
  class CagA       (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp (ref="never or occasionally")
        dm         (ref="no")
        / param=ref;

  model pt_nafld*status(0) =
        CagA
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / eventcode=1 ties=efron;

  hazardratio CagA / diff=ref cl=wald;
run;

proc phreg data=cr_dat;
  class hpcag      (ref='0')
        sex        (ref="female")
        ethnic     (ref="whites")
        smk        (ref="never smoker")
        mvpagrp    (ref="<500 METs")
        drkdosegrp (ref="never or occasionally")
        dm         (ref="no")
        / param=ref;

  model pt_nafld*status(0) =
        hpcag
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / eventcode=1 ties=efron;

  hazardratio hpcag / diff=ref cl=wald;
run;



/**************************************************************************/
/********Covariate interaction analysis: C2-not done         **************/
/**************************************************************************/ 
proc means data=tablecomp n nmiss median min max mean std ; var age tdi bmi; run;

proc freq data=tablecomp; table sex ethnic smk mvpagrp drkdosegrp dm/norow nocol nopercent; run;

data intsample; set tablecomp; 
if age <57 then agecat=0; else if age>=57 then agecat=1;
if tdi <-2.3048100 then tdicat=0; else if tdi >=-2.3048100 then tdicat=1;
if bmi <25 then bmiobs=0; else if bmi>=25 then bmiobs=1;
if smk=0 then eversmk=0; else if smk in (1,2) then eversmk=1;
if mvpagrp in (0,1) then paactive=0; else if mvpagrp=2 then paactive=1; 
if drkdosegrp=0 then everdrk=0; else if drkdosegrp in (1,2,3) then everdrk=1;
run;
proc freq data=intsample; table agecat tdicat bmiobs eversmk paactive everdrk dm/norow nocol nopercent; run;



/**************************************************************************/
/****Firth penalized (bias-reduced) Cox model analysis: B5   **************/
/**************************************************************************/ 
%macro firthcoxmodels(datain=, exposure=, outcome=, time=pt_nafld);
ods output HazardRatios=hr3;
proc phreg data=&datain;
  class &exposure (ref='0')
        sex      (ref="female")
        ethnic   (ref="whites")
        smk      (ref="never smoker")
        mvpagrp  (ref="<500 METs")
		drkdosegrp(ref="never or occasionally")
        dm       (ref="no")
        / param=ref;
  model &time*&outcome(0) =
        &exposure
        age sex ethnic tdi smk mvpagrp drkdosegrp dm bmi
        / ties=BRESLOW firth risklimits=pl;
  hazardratio &exposure / diff=ref cl=pl;
run;
ods output close;
data hr3fmt; set hr3;
Exposure=Description;
HRCI = cat(round(HazardRatio,0.01),' (',round(PLLower,0.01),'-', round(PLUpper,0.01),')');
keep Exposure HRCI;
run;
proc print data=hr3fmt; run;
%mend;

proc freq data=surv; table (hpylori1 CagA hpcag)*nafld/nocol nopercent; run;

%firthcoxmodels(datain=surv,exposure=hpylori1,outcome=nafld); /*N=4246*/
%firthcoxmodels(datain=surv,exposure=CagA,outcome=nafld);
%firthcoxmodels(datain=surv,exposure=hpcag,outcome=nafld);


/**************************************************************************/
/****    Definition of H. pylori seropositivity : B7         **************/
/**************************************************************************/

data diffexp; set surv;
hpcagnew=hpcag;
if hpcag=. then hpcagnew=4;
testhp1=(sumhpy1>=1);
testhp2=(sumhpy1>=2);
testhp3=(sumhpy1>=3);
testhp4=(sumhpy1>=4);
run;

proc freq data=diffexp; table (hpylori1 hpcagnew sumhpy1 testhp:)*nafld/norow nopercent; run;

%coxmodels(datain=diffexp,exposure=hpcagnew,outcome=nafld);

%coxmodels(datain=diffexp,exposure=testhp1,outcome=nafld);
%coxmodels(datain=diffexp,exposure=testhp2,outcome=nafld);
%coxmodels(datain=diffexp,exposure=testhp3,outcome=nafld);
%coxmodels(datain=diffexp,exposure=testhp4,outcome=nafld);


################################################################################
###                  Zhao LG, Mar 1, 2026 CEBP R1                            ###
###                  Hp and NAFLD--meta analysis                             ###
################################################################################

################################################################################
#1 Preparation for database#
################################################################################

################################################################################
##1.1 libraries and self-made functions required
################################################################################
packagelist<-c("reshape","openxlsx","dosresmeta","rms","metafor","sjPlot",
               "sjmisc","readxl","dplyr","forestplot")
sapply(packagelist,library,character.only=T)

getwd()
setwd("C:/Users/lz592/OneDrive - Yale University/NHANES/HpSLD/MS meta/Hp revise_20251121/Submission/5CEBP/CEBP R1_20260209/AnalysisNew")
################################################################################
##1.2 set database, load data, variable calculated
################################################################################
getwd()
workbook <-"Results_MASLD_20260301.xlsx"

alldata <- read_excel(workbook, sheet = "F2dt")
#data process
alldata$logrr <- log(alldata$RR)
alldata$loglb <- log(alldata$LL)
alldata$logub <- log(alldata$UL)
alldata$se <- ((alldata$logub-alldata$loglb)/(2*qnorm(0.975,0,1)))
alldatause <- alldata

################################################################################
##2.1 pooled analyses
################################################################################

#2.1.1 Pooled data
#pdf("F1_20260301.pdf") 
#2.1.1.1 all other studies
alldatauseless <- subset(alldatause, Country != "United Kingdom")
res <- rma(yi = logrr, sei = se, data = alldatauseless, method = "DL")
predict(res,transf=exp)
forest(res, transf=exp, refline=1, mlab="", addcred=TRUE, 
       slab=paste(alldatauseless$Author, alldatauseless$Year, sep=", "))
text(-0.5, -1, pos=4, bquote(
  paste("All Studies (",I^2,"=",.(formatC(res$I2, digits=1, format="f")),"%)")))
text(-0.5,9.5, "Author, Year",pos=4)
text(2.08,9.5, "RR [95% CI]",pos=4)
title("Hp and NAFLD in prospective cohorts",font.main=1,cex.main=1)
#dev.off()

#2.1.1.1 all studies included
res <- rma(yi = logrr, sei = se, data = alldatause, method = "DL")
predict(res,transf=exp)
forest(res, transf=exp, refline=1, mlab="", addcred=TRUE, 
       slab=paste(alldatause$Author, alldatause$Year, sep=", "))
text(-0.5, -1, pos=4, bquote(
 paste("All Studies (",I^2,"=",.(formatC(res$I2, digits=1, format="f")),"%)")))
text(-0.5,9.5, "Author, Year",pos=4)
text(2.08,9.5, "RR [95% CI]",pos=4)
title("Hp and NAFLD in prospective cohorts",font.main=1,cex.main=1)
#dev.off()

#subgroup
alldatause$Serology <- factor(alldatause$Serology)

for (g in levels(alldatause$Serology)) {
  dat <- subset(alldatause, Serology == g)
  res <- rma(yi = logrr, sei = se, data = dat, method = "DL")
  cat("\nSerology:", g, "\n")
  print(predict(res, transf = exp))
}

#2.1.2 publication bias
# Funnel plot
funnel(
  res,
  xlab = "Log Relative Risk",
  ylab = "Standard Error",
  main = "Funnel Plot of Studies on H. pylori and MASLD",
  refline = 0
)

# Egger’s regression test for small-study effects
egger_test <- regtest(res, model = "rma", predictor = "sei")
egger_test

#2.1.3 sensitivity analysis
# Perform leave-one-out analysis
library(metafor)

# Run the analysis
loo_res <- leave1out(res)
loo_df  <- as.data.frame(loo_res)

# Convert to RR scale
loo_df$RR <- exp(loo_df$estimate)
loo_df$LL <- exp(loo_df$ci.lb)
loo_df$UL <- exp(loo_df$ci.ub)

# Order by RR for a tidy plot
ord <- order(loo_df$RR)
loo_df <- loo_df[ord, ]

# Plot
par(mar = c(5, 10, 3, 2))
plot(loo_df$RR, seq_along(loo_df$RR),
     xlim = c(1.1, 1.4),
     xlab = "Pooled Relative Risk (leave-one-out) with 95% CI",
     ylab = "", yaxt = "n", pch = 19)
segments(loo_df$LL, seq_along(loo_df$RR), loo_df$UL, seq_along(loo_df$RR))
abline(v = exp(coef(res)), lty = 2)  # overall pooled RR
axis(2, at = seq_along(loo_df$RR),
     labels = paste("Omitting", paste(alldatause$Author, alldatause$Year, sep = ", "))[ord],
     las = 1, cex.axis = 0.8)
title("Leave-One-Out Sensitivity Analysis")

################################################################################
##2.2 Forest plot for total analysis
################################################################################
 
usetest <- read_excel(workbook, sheet = "F2metdt")

usetest <- usetest %>%
  mutate(Author_Year = ifelse(is.na(Year), "Pooled", paste0(Author, " (", Year, ")")),
         Participants = formatC(Participants, format = "d", big.mark = ","))

tabletext <- rbind(
  c("Author (Year)", "Country", "Cases", "Participants", "RR (95% CI)"),
  cbind(usetest$Author_Year, usetest$Country, usetest$Cases, usetest$Participants, usetest$RRCI)
)

tabletext[9, 1]  <- "Pooled (excluding Zhao, 2025)"
tabletext[11, 1] <- "Overall pooled"

finalplot <- forestplot(
  tabletext,
  hrzl_lines = list("2" = gpar(lty=1, col="gray60", lwd=1)),
  mean  = c(NA, usetest$RR),
  lower = c(NA, usetest$LL),
  upper = c(NA, usetest$UL),
  graph.pos = 5,
  graphwidth = unit(50, "mm"),
  colgap = unit(1, "mm"),
  boxsize = c(NA, rep(0.2, 7), 0.2, 0.2, 0.2),
  clip = c(0.5, 2.5),
  col = fpColors(box="gray36", lines="gray36"),
  vertices = FALSE,
  zero = 1, xticks = c(0.8, 1.0, 1.5, 2.0),
  txt_gp = fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.9), cex=1),
  is.summary = { x <- rep(FALSE, 11); x[c(9, 11)] <- TRUE; x }
)
finalplot

pdf("F1_20260302.pdf", width = 8, height = 4) 
finalplot
dev.off()


