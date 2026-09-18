package stelarx.tree;

/** Newick-safe formatting for logical taxon names. */
public final class NewickLabel {
    private NewickLabel() {}

    /**
     * Keep ordinary identifier labels (including underscores) unquoted. Quote
     * labels containing whitespace or Newick syntax, doubling embedded single
     * quotes according to the Newick convention.
     */
    public static String encode(String label) {
        if (label == null) throw new IllegalArgumentException("Taxon label is null");
        boolean quote = label.isEmpty();
        for (int i = 0; i < label.length() && !quote; i++) {
            char c = label.charAt(i);
            quote = Character.isWhitespace(c)
                || c == '(' || c == ')' || c == '[' || c == ']'
                || c == ',' || c == ':' || c == ';'
                || c == '\'' || c == '"';
        }
        if (!quote) return label;
        return "'" + label.replace("'", "''") + "'";
    }
}
