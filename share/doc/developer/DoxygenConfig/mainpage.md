CFS++                                                              {#mainpage}
=====

Welcome to the  documentation of the <em>Coupled Field System  in C++</em>, or
short <em>CFS++</em>!

Description
-----------

Manuals
-------

### [Developer's Manual](devmanual.html)
### [Markdown Syntax Manual](mdsyntax.html)
### [Coding Rules](codingrules.html)
### [XML Reference](xmlref.html)
### [Schema Documentation](schemadoc.html)

Acknowledgements
----------------

Tips for Conversion from LaTeX to Markdown
------------------------------------------

### Conversion of `.tex` files

LaTeX    can    be   converted    to    Markdown    syntax   using    `pandoc`
(http://johnmacfarlane.net/pandoc). It  is best to convert  a document chapter
by chapter instead of converting the main document, which includes the chapter
`.tex` files. The output of `pandoc` is  assumed to be UTF-8 encoded for further
processing by doxygen.

\warning Note, that LaTeX footnotes, equations and references to them, figures
and code listings are not handled by Doxygen automatically and must be handled
manually!

To control how LaTeX sections get converted to Markdown sections the 
`--base-header-level=N`

In order to get citation references for  Doxygen right, one has to add a extra
`sed` command.

Use

    pandoc -f latex -t markdown -s in.tex | sed -e 's/@/\\cite /g' \ 
                                                -e 's/\\$//'
                                                -e 's/\[\^\([0-9]\)]/\\latexonly\\footnotemark\[\1\]\\endlatexonly\\htmlonly<sup id="footnotemark_\1"><a href="#footnote_\1">\1<\/a><\/sup>\\endhtmlonly/g' | \
                                            iconv -f UTF-8 -t ASCII --unicode-subst='_UNICODE_' | \                                             iconv -f ASCII -t UTF-8 > out.md 

Footnotes

\code
    \latexonly\footnotetext[1]{Agile Software Development: \url{http://www.agile-process.org}}\endlatexonly
    \htmlonly<sup id="footnote_1">1</sup> Agile Software Development: <a href="#footnotemark_1">back up &uarr;</a>\endhtmlonly
\endcode

Images

\code
Introduction{#georg}
------------


\image html images/cfs-logo.png
\htmlonly
<center id="cfs-logo">Figure: Das ist das CFS++ Logo.</center>
\endhtmlonly

\latexonly
\begin{DoxyImage}
\includegraphics[width=7cm]{images/cfs-logo.png}
\caption{Das ist das C\-F\-S++ Logo.}
\label{cfs-logo}
\end{DoxyImage}
\endlatexonly
\endcode

\page codingrules Coding Rules
\htmlinclude coding_rules.html
