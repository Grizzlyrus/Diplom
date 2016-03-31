package utils;

/**
 * Created by Admin on 08.03.2016.
 */
public class Pair<A extends Comparable<? super A>, B extends Comparable<? super B>> implements Comparable<Pair<A, B>> {
    public A first;
    public B second;

    public Pair(Pair<A, B> pair) {
        this.first = pair.first;
        this.second = pair.second;
    }

    public Pair(A first, B second) {
        this.first = first;
        this.second = second;
    }

    public static <A extends Comparable<? super A>, B extends Comparable<? super B>> Pair<A, B> of(A first, B second) {
        return new Pair(first, second);
    }

    public int compareTo(Pair<A, B> o) {
        int cmp = o == null?1:this.first.compareTo(o.first);
        return cmp == 0?this.second.compareTo(o.second):cmp;
    }

    public int hashCode() {
        return 31 * hashcode(this.first) + hashcode(this.second);
    }

    private static int hashcode(Object o) {
        return o == null?0:o.hashCode();
    }

    public boolean equals(Object obj) {
        return !(obj instanceof Pair)?false:(this == obj?true:this.equal(this.first, ((Pair)obj).first) && this.equal(this.second, ((Pair)obj).second));
    }

    private boolean equal(Object o1, Object o2) {
        return o1 == o2 || o1 != null && o1.equals(o2);
    }

    public String toString() {
        return "(" + this.first + ", " + this.second + ')';
    }

    public Pair<A, B> clone() {
        return new Pair(this);
    }
}

