import { useState } from 'react'
import { Button } from '@/components/ui/button'

function App() {
  const [count, setCount] = useState(0)

  return (
    <div className="flex min-h-screen items-center justify-center bg-background">
      <div className="space-y-6 text-center">
        <h1 className="text-4xl font-bold tracking-tight">
          Tailwind CSS + shadcn/ui
        </h1>
        <p className="text-muted-foreground">
          Your setup is complete and working! 🎉
        </p>
        <div className="flex gap-4 justify-center">
          <Button onClick={() => setCount((count) => count + 1)}>
            Count is {count}
          </Button>
          <Button variant="outline">Outline Button</Button>
          <Button variant="secondary">Secondary</Button>
        </div>
        <p className="text-sm text-muted-foreground">
          Edit <code className="bg-muted px-2 py-1 rounded">src/App.tsx</code> to get started
        </p>
      </div>
    </div>
  )
}

export default App
