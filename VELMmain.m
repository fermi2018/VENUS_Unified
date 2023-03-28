 clc
clear
% close all
colordef black
dbstop if error


	addpath('generageomBar','-end')
    addpath('OtticoBar','-end')
    addpath('OtticoBar\new17Optica','-end')
	
load VELM
	
% fil_str='C:\Users\albig\Politecnico\Dottorato\3a_NEGF_VCSEL\Codici\Matlab\NEGF-DD\202008023-D1ANA_TJ_Tibaldi\20220513_VENUS-etchedTJ\dati\MarkusN_BTJvera.str';

                if isfield(mode,'quasi1D') && mode.quasi1D==1
                    [velm] = f_CallVELM_1D(mesh,mode,mode1,ParVet,VelmOptions,fil_str,mode.quasi1D);
                else
                    [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
                end
	