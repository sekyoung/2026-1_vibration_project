%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% University/School Laboratory Report
% LaTeX Class
% Version 4.0 (March 21, 2022)
%
% This class originates from:
% https://www.LaTeXTemplates.com
%
% Authors:
% Vel (vel@latextemplates.com)
% Linux and Unix Users Group at Virginia Tech Wiki 
%
% License:
% CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------------------------------
%	CLASS CONFIGURATION
%----------------------------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{CSUniSchoolLabReport}[2022/03/21 Creodocs University/School Laboratory Report Class v4.0]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}} % Pass through any extra options specified to the base class
\ProcessOptions\relax % Process class options

\LoadClass{article} % Load the base class

%----------------------------------------------------------------------------------------
%	REQUIRED PACKAGES AND MISC CONFIGURATIONS
%----------------------------------------------------------------------------------------

\usepackage{enumitem} % Required for list customization
\setlist{partopsep=0pt, topsep=0pt} % Customize spacing around and inside lists
\renewcommand{\labelenumi}{\alph{enumi}.} % Change numbering in the enumerate environment by letter rather than number

\setlength{\parindent}{0pt} % Suppress paragraph indentation

\usepackage{graphicx} % Required for including images
\graphicspath{{Figures/}{./}} % Specifies where to look for included images (trailing slash required)

\usepackage{float} % Allows more precisely positioning floats e.g. \begin{figure}[H]

\usepackage{mhchem} % Package for chemical equation typesetting
\usepackage{siunitx} % Provides the \SI{}{} and \si{} commands for typesetting technical/scientific SI units correctly

\usepackage{amsmath} % Required for some math elements 
\usepackage{tablefootnote} % Required for footnotes in tables
\usepackage{kotex} % Required for Korean characters
\usepackage{caption} % Add this to your preamble
\usepackage{booktabs} % Required for better horizontal rules in tables
\usepackage{textcomp} % Provides extra symbols, such as \textdegree for the degree symbol
\usepackage{hyperref} % Required for hyperlinks in the PDF
\usepackage{cleveref} % Required for clever referencing of figures, tables, etc.
\usepackage{subcaption} % Required for creating subfigures
\usepackage{pdflscape} % Rotates the page horizontally in the PDF reader
\usepackage{multicol}   % Allows flexible multi-column formatting
%----------------------------------------------------------------------------------------
%	MARGINS
%----------------------------------------------------------------------------------------

\usepackage[
	top=2.54cm, % Top margin
	bottom=2.54cm, % Bottom margin
	left=1.91cm, % Left margin
	right=1.91cm, % Right margin
	%showframe % Uncomment to show frames around the margins for debugging purposes
]{geometry}

%----------------------------------------------------------------------------------------
%	FONTS
%----------------------------------------------------------------------------------------

\usepackage[T1]{fontenc} % Output font encoding for international characters
\usepackage[utf8]{inputenc} % Required for inputting international characters
\usepackage[scale=0.96]{XCharter} % Use the XCharter text font
\usepackage[xcharter]{newtxmath} % Set the math font

%----------------------------------------------------------------------------------------
%	BIBLIOGRAPHY
%----------------------------------------------------------------------------------------

\usepackage{biblatex}
