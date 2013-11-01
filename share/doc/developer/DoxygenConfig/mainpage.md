CFS++                                                              {#mainpage}
=====

Welcome to the  documentation of the <em>Coupled Field System  in C++</em>, or
short <em>CFS++</em>!

Description
-----------

### Umstellung auf Doxygen als einziges Dokutool

- Vorteile
    + Ein  einziges Tool  mit einheitlicher  Syntax sowohl  für Code  als auch
      Manuals	
    + Markdown Syntax für lange Dokumente wird unterstützt
    +      Formeln      in      LaTeX      syntax      (in      HTML      über
      [MathJax](http://www.mathjax.org/))                              möglich
      $\\sum\\limits_{i=1}^{N}\\sqrt{1/i}$	
    + Literaturverzeichnisse über BibTex möglich
    + Ausgabe sowohl in HTML als auch LaTeX/PDF
    + Scrollbare Graphviz Vererbungsdiagramme in SVG

- Nachteile
    + Verweise  (z.B. `\ref` $\\rightarrow$ `\label`)  innerhalb des Dokuments
      schwierig	
    + Code Listings mit Syntax Highlighting nur für wenige Sprachen

- Lösungsmöglichkeiten
    + Verweise auf Bilder, Formeln, Listings und Fussnoten über `INPUT_FILTER`
      bzw. `FILTER_PATTERNS`  und CMake `CONFIGURE_FILE()` in  Kombination mit
      Marken     wie    `@FIG_LABEL_FIGNAME@`,     `@FIG_REF_FIGNAME@`    oder
      `@EQ_LABEL_EQNAME@`, `@EQ_REF_EQNAME@` realisieren.	
    + `\href` in LaTeX Figures auf `\ref` ersetzen.
    + Für HTML MathJax für Formeln lokal installieren	
    + Im  CMake Build Skript  für die Doku  auf Doxygen $\\geq$  1.8.4 testen,
      bzw.  selber  bauen.  Graphviz  `dot`  $\\geq$  2.28  mit  png  und  svg
      Unterstützung.	
    +  Listings über  über  Fenced  Code Blocks  und  dann per  `INPUT_FILTER`
      bzw.     `FILTER_PATTERNS`      durch     Pygments      pipen     (siehe
      [Fenced Code Blocks in der Doxygen Doku](http://www.stack.nl/~dimitri/doxygen/manual/markdown.html#md_fenced)  bzw.  [Pandoc Fenced  Code
      Blocks](http://johnmacfarlane.net/pandoc/README.html#fenced-code-blocks)).
~~~~ {#mycode .c}
#include <stdio.h>

int main(int argc, char** argv) {
  std::cout << "Hello from C++!" << std::endl;
  return 0;
}
~~~~

### Schema Dokumentation
- Schema über  `<annotation>` Tags dokumentieren  hat den Vorteil,  dass einem
  kontextsensitive Hilfe im XML Editor angezeigt wird.
~~~~ {.xml}
<xsd:annotation>
  <xsd:documentation xml:lang="en">
    Coupled Field Solver project CFS++
    Schema for basic PDE description and boundary conditions
    All other PDE descriptions are derived from this type
  </xsd:documentation>
</xsd:annotation>
~~~~
- Schema    Dokumentation    kann    am    einfachsten    und    besten    mit
  [Oxygen nach HTML und PDF generiert](http://www.oxygenxml.com/doc/ug-editor/topics/documentation-XML-Schema.html)
  werden (siehe [Schema Doku auf wiki.mdmt](http://wiki.mdmt.tuwien.ac.at/cfs_docu/schema_sim)).

### Tips for Conversion from LaTeX to Markdown

#### Conversion of `.tex` files

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

#### Footnotes

\code
    \latexonly\footnotetext[1]{Agile Software Development: \url{http://www.agile-process.org}}\endlatexonly
    \htmlonly<sup id="footnote_1">1</sup> Agile Software Development: <a href="#footnotemark_1">back up &uarr;</a>\endhtmlonly
\endcode

#### Images

\code

    # The FILTER_PATTERNS tag can be used to specify filters on a per file pattern
    # basis.  Doxygen will compare the file name with each pattern and apply the
    # filter if there is a match.
    #
    # The filters are a list of the form:
    #   pattern=filter (like *.cpp=my_cpp_filter).
    # See INPUT_FILTER for further info on how filters are used. If
    # FILTER_PATTERNS is empty or if non of the patterns match the file name,
    # INPUT_FILTER is applied.
    FILTER_PATTERNS = *.cpp="awk 'BEGIN {print "#define __cplusplus"} {print}'" \
                      *.hpp="awk 'BEGIN {print "#define __cplusplus"} {print}'"
\endcode


#### [HTML Documentation with Equation Numbers - Referencing an External PDF Document with Doxygen’s HTML Documentation](http://www.cheshirekow.com/wordpress/?p=335)


Manuals
-------

### [Developer's Manual](devmanual.html)
### Markdown
#### [Markdown Syntax Manual](mdsyntax.html)
#### [Markdown Support in Doxygen](http://www.stack.nl/~dimitri/doxygen/manual/markdown.html)
#### Markdown Editors
* [MarkdownPad is a full-featured Markdown editor for Windows.](http://markdownpad.com/)
* [Markdown Mode for Emacs](http://jblevins.org/projects/markdown-mode/)
* [Markdown Plugin](http://www.norberteder.com/blog/post/Markdown-Unterstutzung-fur-Sublime-Text.aspx) [Sublime Text](http://www.sublimetext.com/)
* [Editorial (iPad)](http://omz-software.com/editorial/)
* [iA Writer (Mac/iPad/iPhone)](http://www.iawriter.com/)
* [Lightpaper (Mac/Android)](http://clockworkengine.com/)
* [Draft (Android)](http://www.mvilla.it/draft/)

### [Coding Rules](codingrules.html)
### [XML Reference](xmlref.html)
### [Schema Documentation](schemadoc.html)

Acknowledgements
----------------

