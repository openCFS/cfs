Doxygen Documentation {#mainpage}
=================================

Markdown & Doxygen
------------------

See the [doxygen-manual](https://www.stack.nl/~dimitri/doxygen/manual/markdown.html) for further information.
There are some important points:
* The page should be started with a level 1 header including a [header id](https://www.stack.nl/~dimitri/doxygen/manual/markdown.html#md_header_id) like `# Level 1 Heading {#headerID}`.
* general pages should be placed in `share/doc/developer/doxygen/pages`

To create links between pages use `{#pageid}` in one page.
This defines the name of the created html file.
You can then link to this page by a link like

    [page YXZ](pageid.html)

Here are some exapmles of links within doxygen to
* a [specific page](Markdown_Syntax.html) or its [header id](@ref Markdown_Syntax)
* a [doxygen html file](classCoupledField_1_1CoefFunction.html)
* elements documented by doxygen like [AcousticPDE](@ref CoupledField::AcousticPDE) or [GDataInfo](@ref GDataInfo)
* for CFS-classes do not forget to include the [namespace](@ref CoupledField), i.e. `[linkname](CoupledField::AcousticPDE)`

For further information see [doxygen markdown documentation](http://www.stack.nl/~dimitri/doxygen/manual/markdown.html).

Snippets
----------------
Doxygen allows to reference code snippets in the documentation.
This is acomplised by placing markers like `//! [markertext]` around a piece of code.
The marker must appear *exactly twice*, and it must be placed at the *beginning of a line* (no whitespace).
It may then be referenced by `\snippet filename markertext`.

For the following example the marker `//! [Read PtrParamNode]` has been inserted into the file `source/DataInOut/ParamHandling/XMLMaterialHandler.cc`.
The snippet may then be referenced by using `\snippet source/DataInOut/ParamHandling/XMLMaterialHandler.cc Read PtrParamNode` in any documentation file (e.g. this markdown file) resulting in
\snippet XMLMaterialHandler.cc Read PtrParamNode
i.e. the code section is typset in doxygen.
The path to the source file can start from `source`, but does not have to.

For further details see the [sippets section of the doxygen documentation](http://www.stack.nl/~dimitri/doxygen/manual/commands.html#cmdsnippet).
